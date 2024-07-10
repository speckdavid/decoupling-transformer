#include "mwis_factoring.h"

#include "../algorithms/max_cliques.h"
#include "../algorithms/sccs.h"

#include "../plugins/plugin.h"

#include "../utils/timer.h"
#include "../utils/hash.h"

#include "../task_utils/causal_graph.h"
#include "../task_proxy.h"
#include "../tasks/root_task.h"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <queue>


using namespace std;


namespace decoupling {
void MWISFactoring::PotentialLeaf::add_leaf_only_schema(size_t as_id, const ActionSchema &action_schema) {
    if (std::find(action_schemes.begin(),
                  action_schemes.end(),
                  as_id) == action_schemes.end()){
        action_schemes.push_back(as_id);
        num_actions += action_schema.num_actions;
        bool all_pre_vars_in_vars = true;
        for (int pre_var : action_schema.pre_vars){
            if (!std::binary_search(vars.begin(), vars.end(), pre_var)){
                all_pre_vars_in_vars = false;
                break;
            }
        }
        if (all_pre_vars_in_vars){
            self_mobile_as.push_back(as_id);
        }
    }
}

MWISFactoring::MWISFactoring(const plugins::Options &opts) : Factoring(opts),
        strategy(opts.get<WMIS_STRATEGY>("strategy")),
        min_mobility(opts.get<int>("min_mobility")),
        min_flexibility(opts.get<double>("min_flexibility")),
        min_fact_flexibility(opts.get<double>("min_fact_flexibility")),
        add_cg_sccs_(opts.get<bool>("add_cg_sccs")) {

    if (log.is_at_least_normal()) {
        log << "Using LP factoring with strategy: ";
        switch (strategy) {
            case WMIS_STRATEGY::MML:
                log << "maximize number of mobile leaves." << endl; break;
            case WMIS_STRATEGY::MMAS:
                log << "maximize number of mobile action schemas." << endl; break;
            case WMIS_STRATEGY::MM_OPT:
                log << "maximize leaf mobility (exact)." << endl; break;
            case WMIS_STRATEGY::MFA:
                log << "maximize number of mobile facts." << endl; break;
            case WMIS_STRATEGY::MM:
                log << "maximize leaf mobility (sum)." << endl; break;
            case WMIS_STRATEGY::MCL:
                log << "maximize number of mobile conclusive leaves." << endl; break;
            case WMIS_STRATEGY::MCM:
                log << "maximize conclusive leaf mobility." << endl; break;
            default:
                log << "ERROR: unknown LP factoring strategy." << endl;
                exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
        }
    }

    if (min_number_leaves > 1 && strategy != WMIS_STRATEGY::MML){
        log << "WARNING: WMIS factoring does not support setting a minimal number of leaf factors." << endl;
        min_number_leaves = 1;
    }

    if (strategy != WMIS_STRATEGY::MFA && min_fact_flexibility > 0.0){
        log << "Option min_fact_flexibility is only possible in combination with strategy MFA." << endl;
        exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }
}

inline double get_log(double num_actions) {
    return log(max(1.0001, num_actions));
}

void MWISFactoring::add_leaf_intersection_edges(Graph &graph,
                                                const vector<vector<size_t>> &var_to_p_leaves) const {
    // non-empty intersection between potential leaves
    vector<vector<bool> > pleaf_intersect(potential_leaf_nodes.size() - 1);
    for (size_t i = 0; i < pleaf_intersect.size(); ++i){
        pleaf_intersect[i].resize(potential_leaf_nodes.size() - i - 1, false);
    }
    for (int var = 0; var < (int) task->get_num_variables(); ++var){
        for (size_t pot_leaf_1 : var_to_p_leaves[var]){
            for (size_t pot_leaf_2 : var_to_p_leaves[var]){
                if (pot_leaf_1 < pot_leaf_2){
                    pleaf_intersect[pot_leaf_1][pot_leaf_2 - pot_leaf_1 - 1] = true;
                }
            }
        }
    }
    for (size_t p_leaf_1 = 0; p_leaf_1 < potential_leaf_nodes.size(); ++p_leaf_1) {
        for (size_t p_leaf_2 = p_leaf_1 + 1; p_leaf_2 < potential_leaf_nodes.size(); ++p_leaf_2) {
            // we need p_leaf_1 < p_leaf_2
            if (pleaf_intersect[p_leaf_1][p_leaf_2 - p_leaf_1 - 1]) {
                graph[p_leaf_1].push_back(p_leaf_2);
                graph[p_leaf_2].push_back(p_leaf_1);
            }
        }
        // TODO: is this actually needed?
        utils::sort_unique(graph[p_leaf_1]);
    }
}

void MWISFactoring::add_outside_pre_var_edges(Graph &graph,
                                              const vector<vector<size_t>> &var_to_p_leaves) const {
    for (size_t i = 0; i < potential_leaf_nodes.size(); ++i){
        const PotentialLeafNode &pleaf = potential_leaf_nodes[i];
        for (int var : pleaf.outside_pre_vars){
            for (size_t pleaf_id : var_to_p_leaves[var]){
                graph[i].push_back(pleaf_id);
                graph[pleaf_id].push_back(i);
            }
        }
    }
    for (size_t i = 0; i < potential_leaf_nodes.size(); ++i){
        // TODO avoid this
        utils::sort_unique(graph[i]);
    }
}

bool MWISFactoring::is_as_leaf_irrelevant(const ActionSchema &as, const PotentialLeaf &leaf) const {
    assert(!has_as_pre_or_eff_on_leaf(as, leaf));
    vector<bool> is_leaf_pre_var(task->get_num_variables(), false);
    for (auto as_id : leaf.action_schemes){
        for (int var : action_schemas[as_id].pre_vars){
            is_leaf_pre_var[var] = true;
        }
    }
    for (int evar : as.eff_vars){
        if (is_leaf_pre_var[evar]){
            return false;
        }
    }
    return true;
}

bool MWISFactoring::is_as_leaf_conclusive(const ActionSchema &as, const PotentialLeaf &leaf) {
    // is every var in leaf contained in either the precondition of effect of as
    for (int lvar : leaf.vars) {
        bool covered = false;
        for (int pvar: as.pre_vars) {
            if (lvar == pvar){
                covered = true;
            } else if (pvar > lvar){
                break;
            }
        }
        if (!covered) {
            for (int evar: as.eff_vars) {
                if (lvar == evar) {
                    covered = true;
                    break;
                } else if (evar > lvar) {
                    return false;
                }
            }
            if (!covered){
                return false;
            }
        }
    }
    return true;
}

bool MWISFactoring::has_as_pre_or_eff_on_leaf(const ActionSchema &as, const PotentialLeaf &leaf) {
    for (int lvar : leaf.vars) {
        for (int pvar: as.pre_vars) {
            if (lvar == pvar){
                return true;
            } else if (pvar > lvar){
                break;
            }
        }
        for (int evar: as.eff_vars) {
            if (lvar == evar){
                return true;
            } else if (evar > lvar){
                break;
            }
        }
    }
    return false;
}

void MWISFactoring::construct_graph_conclusive_leaves(Graph &/*graph*/) {

    // TODO implement this
    cerr << "not implemented in mwis_factoring.cc" << endl;
    utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);

//    assert(variables.size() == 0);
//    assert(constraints.size() == 0);
//
//    compute_action_schemas();
//
//    if (action_schemas.empty()){
//        // mostly for trivially unsolvable task from translator?
//        log << "ERROR: No action schemas." << endl;
//        return;
//    }
//
//    compute_potential_leaves();
//
//    if (static_cast<int>(potential_leaves.size()) < min_number_leaves){
//        log << "Only " << potential_leaves.size() <<
//            " potential leaves left, but minimum number of leaves is " <<
//            min_number_leaves << "." << endl;
//        return;
//    }
//
//    if (!check_timeout()){
//        return;
//    }
//
//    vector<bool> can_be_leaf_var(task->get_num_variables(), false);
//    // binary var for each potential leaf; 1 => becomes leaf, 0 => not
//    for (const auto &pleaf : potential_leaves) {
//        for (int var : pleaf.vars){
//            can_be_leaf_var[var] = true;
//        }
//        double obj = 0.0;
//        if (strategy == WMIS_STRATEGY::MCL) {
//            // for tie-breaking, give objective value of 1 to every mobile leaf
//            obj = 1.0;
//        }
//        variables.push_back(lp::LPVariable(0.0, 1.0, obj, true));
//    }
//
//    // an LP variable for every potential leaf that is 1 iff the leaf is conclusive
//    vector<int> concl_pleaf_lp_vars;
//    if (strategy == WMIS_STRATEGY::MCL) {
//        assert(variables.size() == static_cast<int>(potential_leaves.size()));
//        for (size_t i = 0; i < potential_leaves.size(); ++i) {
//            concl_pleaf_lp_vars.push_back(variables.size());
//            variables.push_back(lp::LPVariable(0.0, 1.0, 100.0, true));
//        }
//        assert(concl_pleaf_lp_vars.size() == potential_leaves.size());
//        assert(variables.size() == 2 * static_cast<int>(potential_leaves.size()));
//    }
//
//    // binary var for each action schema of each potential leaf; 1 => is mobile, 0 => not
//    vector<vector<size_t>> mob_as_var_ids(potential_leaves.size());
//    // binary var for each action schema of each potential leaf; 1 => is conclusive for potential leaf, 0 => not
//    vector<vector<size_t>> concl_as_var_ids(potential_leaves.size());
//
//    // for every potential leaf we store the list of action schemas relevant for the constraints;
//    // for MCL, this is the schemas that can be neither irrelevant nor conclusive for the leaf,
//    // for MCM, it's the schemas that become irrelevant or conclusive as global actions
//    vector<vector<size_t>> relevant_as_by_pleaf(potential_leaves.size());
//
//    // list of LP variable IDs that represent mob_as_var_ids for every action schemas, i.e.,
//    // the corresponding LP variable IDs across potential leaves
//    vector<vector<size_t>> as_to_mob_as_var_ids(action_schemas.size());
//
//    for (const PotentialLeaf &pleaf : potential_leaves){
//        for (auto as_id : pleaf.action_schemes){
//            as_to_mob_as_var_ids[as_id].push_back(variables.size());
//            mob_as_var_ids[pleaf.id].push_back(variables.size());
//            variables.push_back(lp::LPVariable(0.0, 1.0, 0.0, true));
//        }
//        // for all action schemas as_id that are not part of pleaf, add a variable that represents if as_id is conclusive for pleaf
//        for (size_t as_id = 0; as_id < action_schemas.size(); ++as_id) {
//            if (find(pleaf.action_schemes.begin(), pleaf.action_schemes.end(), as_id) != pleaf.action_schemes.end()){
//                // action schema is in pleaf => cannot be conclusive for pleaf
//                continue;
//            }
//            bool conclusive_or_irrelevant = false;
//            if (has_as_pre_or_eff_on_leaf(action_schemas[as_id], pleaf)){
//                if (is_as_leaf_conclusive(action_schemas[as_id], pleaf)){
//                    conclusive_or_irrelevant = true;
//                }
//            } else if (is_as_leaf_irrelevant(action_schemas[as_id], pleaf)){
//                conclusive_or_irrelevant = true;
//            }
//            if (strategy == WMIS_STRATEGY::MCM){
//                if (!conclusive_or_irrelevant){
//                    // the action schema cannot be conclusive for this potential leaf, no need to add variables for the objective value
//                    continue;
//                }
//            } else {
//                assert(strategy == WMIS_STRATEGY::MCL);
//                if (conclusive_or_irrelevant){
//                    // the action schema will be conclusive or irrelevant for the leaf, no need to add a constraint
//                    continue;
//                }
//            }
//
//            relevant_as_by_pleaf[pleaf.id].push_back(as_id);
//
//            if (strategy == WMIS_STRATEGY::MCM) {
//                concl_as_var_ids[pleaf.id].push_back(variables.size());
//                variables.push_back(lp::LPVariable(0.0,
//                                                   1.0,
//                                                   action_schemas[as_id].num_actions,
//                                                   true));
//            }
//        }
//    }
//
//    if (min_flexibility > 0.0){
//        add_min_flexibility_constraints(constraints, mob_as_var_ids);
//    }
//
//    if (min_mobility > 1){
//        add_min_mobility_constraints(constraints, mob_as_var_ids);
//    }
//
//    vector<size_t> c_vars_ids = add_center_variables_and_get_ids(variables, constraints, can_be_leaf_var);
//
//    add_leaf_intersection_edges(constraints);
//
//    if (!check_timeout()){
//        return;
//    }
//
//    add_potential_leaf_to_action_schema_constraints(constraints, mob_as_var_ids, can_be_leaf_var, c_vars_ids);
//
//    // add constraints such that conclusiveness of action schemas for pleaves is set properly
//    for (const auto &pleaf : potential_leaves) {
//        if (strategy == WMIS_STRATEGY::MCL) {
//            lp::LPConstraint constraint1(-infty, 0.0);
//            constraint1.insert(pleaf.id, -1.0);
//            constraint1.insert(concl_pleaf_lp_vars[pleaf.id], 1.0);
//            // pleaf can only be conclusive if it becomes a (mobile) leaf
//            constraints.push_back(constraint1);
//        }
//        for (size_t i = 0; i < relevant_as_by_pleaf[pleaf.id].size(); ++i) {
//            size_t as_id = relevant_as_by_pleaf[pleaf.id][i];
//            if (strategy == WMIS_STRATEGY::MCM){
//                lp::LPConstraint constraint1(-infty, 0.0);
//                constraint1.insert(pleaf.id, -1.0);
//                constraint1.insert(concl_as_var_ids[pleaf.id][i], 1.0);
//                // action schema can only be conclusive for pleaf if pleaf becomes a leaf
//                constraints.push_back(constraint1);
//
//                // action schema can be conclusive for pleaf only if it is a global action schema, i.e., iff it is not mobile for any potential leaf
//                for (auto as_mob_var_id : as_to_mob_as_var_ids[as_id]) {
//                    lp::LPConstraint constraint2(-infty, 1.0);
//                    constraint2.insert(concl_as_var_ids[pleaf.id][i], 1.0);
//                    constraint2.insert(as_mob_var_id, 1.0);
//                    constraints.push_back(constraint2);
//                }
//            } else {
//                assert(strategy == WMIS_STRATEGY::MCL);
//
//                // action schema cannot be conclusive or irrelevant for pleaf so needs to be leaf-only for some potential leaf
//                lp::LPConstraint constraint2(-infty, 0.0);
//                constraint2.insert(concl_pleaf_lp_vars[pleaf.id], 1.0);
//                for (auto as_mob_var_id : as_to_mob_as_var_ids[as_id]) {
//                    constraint2.insert(as_mob_var_id, -1.0);
//                }
//                constraints.push_back(constraint2);
//            }
//        }
//    }
}

void MWISFactoring::construct_graph(Graph &graph) {
    assert(graph.empty());

    compute_action_schemas();

    if (action_schemas.empty()){
        // mostly for trivially unsolvable task from translator?
        log << "ERROR: No action schemas." << endl;
        return;
    }

    compute_potential_leaves();

    if (!check_timeout()){
        return;
    }

    if (static_cast<int>(potential_leaf_nodes.size()) < min_number_leaves){
        log << "Only " << potential_leaf_nodes.size() <<
            " potential leaves left, but minimum number of leaves is " <<
            min_number_leaves << "." << endl;
        return;
    }

    graph.resize(potential_leaf_nodes.size());

    vector<vector<size_t>> var_to_p_leaves(task->get_num_variables()); // maps variables to potential leaf ids
    for (size_t i = 0; i < potential_leaf_nodes.size(); ++i) {
        const PotentialLeafNode &pleaf = potential_leaf_nodes[i];
        for (int var : pleaf.vars) {
            var_to_p_leaves[var].push_back(i);
        }
    }

    add_leaf_intersection_edges(graph, var_to_p_leaves);

    if (!check_timeout()){
        return;
    }

    add_outside_pre_var_edges(graph, var_to_p_leaves);
}

vector<int> MWISFactoring::solve_wmis(const Graph &graph,
                                      const vector<double> &weights,
                                      const utils::CountdownTimer &timer) {

    vector<int> independent_set;

    utils::g_log << "Computing max weighted independent set..." << flush;
    double weight = max_cliques::compute_max_weighted_independent_set(graph, weights, independent_set, timer.get_remaining_time());
    utils::g_log << "done!" << endl;

    if (log.is_at_least_verbose()) {
        log << "Weight of computed independent set: " << weight << endl;
    }

    return independent_set;
}

void MWISFactoring::compute_factoring_() {

    // successor node IDs for all graph nodes
    vector<vector<int>> graph;

    if (strategy == WMIS_STRATEGY::MCL || strategy == WMIS_STRATEGY::MCM){
        construct_graph_conclusive_leaves(graph);
    } else {
        construct_graph(graph);
    }

    if (!check_timeout()){
        return;
    }

    if (graph.empty()){
        log << "WARNING: no graph nodes created, stopping." << endl;
        return;
    }

    // save memory
    Factoring::save_memory();

    vector<double> weights(potential_leaf_nodes.size());
    int i = 0;
    for (const auto &pleaf : potential_leaf_nodes){
        weights[i++] = pleaf.weight;
    }

    vector<int> solution = solve_wmis(graph, weights, factoring_timer);

    if (solution.empty()) {
        log << "WARNING: no solution found." << endl;
        return;
    }

    if (static_cast<int>(solution.size()) < min_number_leaves) {
        log << "WARNING: no factoring found with at least " << min_number_leaves << " leaves." << endl;
        return;
    }

    vector<bool> is_leaf_var(task->get_num_variables(), false);
    for (int leaf_id : solution) {
        set<int> leaf;
        for (int var : potential_leaf_nodes[leaf_id].vars) {
            assert(!is_leaf_var[var]);
            is_leaf_var[var] = true;
            leaf.insert(var);
        }
        leaves.emplace_back(leaf.begin(), leaf.end());
    }
}

void MWISFactoring::save_memory() {
    Factoring::save_memory();
    vector<PotentialLeafNode>().swap(potential_leaf_nodes);
    vector<vector<size_t>>().swap(variables_to_action_schemas);
}

inline double compute_leaf_fact_flexibility(const vector<int> &vars,
                                            const AbstractTask &task,
                                            const vector<size_t> &included_as,
                                            const vector<vector<unordered_map<size_t, size_t>>> &facts_to_mobility,
                                            const vector<vector<size_t>> &sum_fact_mobility) {
    double sum_fact_flexibility = 0;
    for (int var: vars) {
        int dom_size = task.get_variable_domain_size(var);
        for (int val = 0; val < dom_size; ++val) {
            double fact_mobility = 0;
            for (size_t as_id: included_as) {
                if (facts_to_mobility[var][val].count(as_id) > 0){
                    fact_mobility += facts_to_mobility[var][val].at(as_id);
                }
            }
            if (sum_fact_mobility[var][val] > 0) {
                // some facts might never be set by any action
                sum_fact_flexibility += fact_mobility / sum_fact_mobility[var][val];
            }
        }
    }
    return sum_fact_flexibility;
}

bool MWISFactoring::fulfills_min_flexibility_and_mobility(
        const PotentialLeaf &pleaf,
        const vector<size_t> &included_as,
        const vector<vector<unordered_map<size_t, size_t>>> &facts_to_mobility,
        const vector<vector<size_t>> &sum_fact_mobility) const {

    if (min_flexibility > 0.0 || min_mobility > 1){
        size_t num_leaf_only_actions = 0;
        for (size_t as : included_as){
            num_leaf_only_actions += action_schemas[as].num_actions;
        }
        if (min_flexibility > 0.0){
            assert(pleaf.num_actions > 0);
            assert(pleaf.num_affecting_actions > 0);
            assert(pleaf.num_affecting_actions >= pleaf.num_actions);
            if (min_flexibility > num_leaf_only_actions / pleaf.num_affecting_actions){
                return false;
            }
        }
        if (min_mobility > 1) {
            if (num_leaf_only_actions > min_mobility){
                return false;
            }
        }
    }
    if (min_fact_flexibility > 0){
        assert(facts_to_mobility.size() == static_cast<size_t>(task->get_num_variables()));
        assert(facts_to_mobility.size() == sum_fact_mobility.size());

        double sum_fact_flexibility = compute_leaf_fact_flexibility(pleaf.vars,
                                                                    *task,
                                                                    included_as,
                                                                    facts_to_mobility,
                                                                    sum_fact_mobility);

        if (min_fact_flexibility > sum_fact_flexibility){
            return false;
        }
    }
    return true;
}

void MWISFactoring::multiply_out_potential_leaf(const vector<pair<vector<int>, vector<size_t>>> &outside_pre_and_ases,
                                                const PotentialLeaf &pleaf,
                                                vector<int> &outside_pre_vars,
                                                vector<size_t> &included_as,
                                                size_t depth,
                                                int &ignored_leaf_candidates,
                                                const vector<vector<unordered_map<size_t, size_t>>> &facts_to_mobility,
                                                const vector<vector<size_t>> &sum_fact_mobility) {
    if (depth == outside_pre_and_ases.size()){
        if (included_as.empty()){
            // no self-mobile AS + no AS included
            return;
        }

        double weight;

        switch (strategy) {
            case WMIS_STRATEGY::MMAS:
                weight = included_as.size();
                break;
            case WMIS_STRATEGY::MFA:
                weight = compute_leaf_fact_flexibility(pleaf.vars,
                                                       *task,
                                                       included_as,
                                                       facts_to_mobility,
                                                       sum_fact_mobility);
                break;
            case WMIS_STRATEGY::MML:
                weight = 1;
                break;
            case WMIS_STRATEGY::MM_OPT:
                weight = 0;
                for (auto as : included_as){
                    weight += action_schemas[as].num_actions;
                }
                weight = get_log(weight);
                break;
            case WMIS_STRATEGY::MM:
                weight = 0;
                for (auto as : included_as){
                    weight += action_schemas[as].num_actions;
                }
                break;
            default: log << "strategy not supported" << endl; exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
        }

        if (fulfills_min_flexibility_and_mobility(pleaf, included_as, facts_to_mobility, sum_fact_mobility)) {
            assert(std::find(potential_leaf_nodes.begin(),
                             potential_leaf_nodes.end(),
                             PotentialLeafNode(outside_pre_vars, pleaf.vars, weight)) == potential_leaf_nodes.end());
            potential_leaf_nodes.emplace_back(outside_pre_vars, pleaf.vars, weight);
        } else {
            ignored_leaf_candidates++;
        }
        return;
    }

    const auto &[outside_pre, ases] = outside_pre_and_ases[depth];
    assert(std::all_of(ases.begin(), ases.end(), [&](size_t as) {
        return find(pleaf.self_mobile_as.begin(), pleaf.self_mobile_as.end(), as) == pleaf.self_mobile_as.end();}));
    size_t num_outside_pre_before = outside_pre_vars.size();
    size_t num_included_as_before = included_as.size();

    // (1) include this AS
    outside_pre_vars.insert(outside_pre_vars.end(), outside_pre.begin(), outside_pre.end());
    included_as.insert(included_as.end(), ases.begin(), ases.end());

    multiply_out_potential_leaf(outside_pre_and_ases, pleaf, outside_pre_vars, included_as, depth + 1, ignored_leaf_candidates, facts_to_mobility, sum_fact_mobility);

    assert(num_included_as_before <= included_as.size());
    included_as.resize(num_included_as_before);
    assert(num_outside_pre_before <= outside_pre_vars.size());
    outside_pre_vars.resize(num_outside_pre_before);

    // (2) do not include this AS
    multiply_out_potential_leaf(outside_pre_and_ases, pleaf, outside_pre_vars, included_as, depth + 1, ignored_leaf_candidates, facts_to_mobility, sum_fact_mobility);
}

void MWISFactoring::multiply_out_action_schemas(
        const vector<PotentialLeaf> &potential_leaves,
        const vector<vector<unordered_map<size_t, size_t>>> &facts_to_mobility,
        const vector<vector<size_t>> &sum_fact_mobility) {
    assert(potential_leaf_nodes.empty());

    if (min_flexibility > 0) {
        compute_variables_to_action_schemas_map();
    }

    int ignored_leaf_candidates = 0;
    for (auto &pleaf: potential_leaves){
        vector<pair<vector<int>, vector<size_t>>> outside_pre_and_as;
        {
            vector<bool> is_leaf_var(task->get_num_variables(), false);
            for (int var: pleaf.vars) {
                is_leaf_var[var] = true;
            }
            utils::HashMap<vector<int>, vector<size_t>> as_by_outside_pre_vars;
            for (size_t as_id: pleaf.action_schemes) {
                if (find(pleaf.self_mobile_as.begin(), pleaf.self_mobile_as.end(), as_id) ==
                    pleaf.self_mobile_as.end()) {
                    vector<bool> is_outside_pre_var(task->get_num_variables(), false);
                    const ActionSchema &as = action_schemas[as_id];
                    for (int var: as.pre_vars) {
                        if (!is_leaf_var[var]) {
                            is_outside_pre_var[var] = true;
                        }
                    }
                    vector<int> outside_pre_vars;
                    for (int var = 0; var < task->get_num_variables(); ++var) {
                        if (is_outside_pre_var[var]) {
                            outside_pre_vars.push_back(var);
                        }
                    }
                    as_by_outside_pre_vars[outside_pre_vars].push_back(as_id);
                }
            }

            outside_pre_and_as.reserve(as_by_outside_pre_vars.size());
            for (const auto &[outside_pre, ases]: as_by_outside_pre_vars) {
                outside_pre_and_as.emplace_back(outside_pre, ases);
            }

            if (min_flexibility > 0){
                size_t num_ops = 0;
                vector<bool> handled_as(action_schemas.size(), false);
                for (int var : pleaf.vars) {
                    for (size_t as : variables_to_action_schemas[var]) {
                        if (!handled_as[as]) {
                            num_ops += action_schemas[as].num_actions;
                            handled_as[as] = true;
                        }
                    }
                }
                pleaf.num_affecting_actions = num_ops;
                if (min_flexibility > pleaf.num_actions / num_ops) {
                    // this is the maximum this leaf can possibly get and it is not enough
                    continue;
                }
            }

            if (min_mobility > pleaf.num_actions) {
                // this is the maximum this leaf can possibly get and it is not enough
                continue;
            }
            if (min_fact_flexibility > 0){
                double max_leaf_fact_flex = compute_leaf_fact_flexibility(pleaf.vars,
                                                                          *task,
                                                                          pleaf.action_schemes,
                                                                          facts_to_mobility,
                                                                          sum_fact_mobility);
                if (min_fact_flexibility > max_leaf_fact_flex){
                    // this is the maximum this leaf can possibly get and it is not enough
                    continue;
                }
            }
        }

        vector<int> outside_pre_vars;
        vector<size_t> included_as(pleaf.self_mobile_as);
        multiply_out_potential_leaf(outside_pre_and_as,
                                    pleaf,
                                    outside_pre_vars,
                                    included_as,
                                    0,
                                    ignored_leaf_candidates,
                                    facts_to_mobility,
                                    sum_fact_mobility);

        if (!check_timeout()){
            return;
        }
    }

    log << "Ignored " << ignored_leaf_candidates << " potential leaves due to minimum mobility / flexibility." << endl;
    log << "Number final leaf candidates: " << potential_leaf_nodes.size() << endl;
}

vector<size_t> MWISFactoring::get_superset_pleaf_ids(const PotentialLeaf &pleaf,
                                                     const vector<PotentialLeaf> &potential_leaves,
                                                     const vector<vector<size_t>> &var_to_p_leaves) {
    vector<size_t> superset_pleaf_ids;
    if (pleaf.vars.size() == 1) {
        for (size_t index : var_to_p_leaves[pleaf.vars[0]]) {
            superset_pleaf_ids.push_back(index);
        }
    } else {
        // since we consider all supersets of p_leaf, we can simply take p_leaf.vars[0],
        // since all supersets must also include that variable
        for (size_t check_pleaf_id : var_to_p_leaves[pleaf.vars[0]]) {
            bool superset_schema = true;
            for (int var : pleaf.vars) {
                if (!binary_search(potential_leaves[check_pleaf_id].vars.begin(), potential_leaves[check_pleaf_id].vars.end(), var)){
                    superset_schema = false;
                    break;
                }
            }
            if (superset_schema) {
                superset_pleaf_ids.push_back(check_pleaf_id);
            }
        }
    }
    return superset_pleaf_ids;
}

void MWISFactoring::compute_variables_to_action_schemas_map() {
    if (variables_to_action_schemas.empty()) {
        variables_to_action_schemas.resize(task->get_num_variables());
        for (size_t as_id = 0; as_id < action_schemas.size(); ++as_id) {
            for (int var: action_schemas[as_id].eff_vars) {
                variables_to_action_schemas[var].push_back(as_id);
            }
        }
    }
}

void MWISFactoring::compute_fact_flexibility(
        vector<vector<unordered_map<size_t, size_t>>> &facts_to_mobility,
        vector<vector<size_t>> &sum_fact_mobility) {
    assert(facts_to_mobility.empty());
    assert(sum_fact_mobility.empty());

    VariablesProxy variables = task_proxy.get_variables();

    // need to store information for effect *facts*, which gets otherwise lost
    facts_to_mobility.resize(task->get_num_variables());
    for (int var = 0; var < task->get_num_variables(); ++var){
        facts_to_mobility[var].resize(variables[var].get_domain_size());
    }

    compute_variables_to_action_schemas_map();

    for (OperatorProxy op : task_proxy.get_operators()) {
        for (EffectProxy eff : op.get_effects()) {
            int eff_var = eff.get_fact().get_variable().get_id();
            for (size_t as_id : variables_to_action_schemas[eff_var]){
                facts_to_mobility[eff_var][eff.get_fact().get_value()][as_id]++;
            }
        }
    }

    sum_fact_mobility.resize(task->get_num_variables());
    for (int var = 0; var < task->get_num_variables(); ++var){
        sum_fact_mobility[var].resize(variables[var].get_domain_size());
        for (int val = 0; val < variables[var].get_domain_size(); ++val){
            for (const auto& [as, num] : facts_to_mobility[var][val]){
                sum_fact_mobility[var][val] += num;
            }
        }
    }
}

void MWISFactoring::compute_potential_leaves() {
    assert(!action_schemas.empty());

    vector<PotentialLeaf> potential_leaves;

    {
        VariablesProxy variables = task_proxy.get_variables();
        utils::HashMap<vector<int>, size_t> leaf_lookup;
        for (size_t as = 0; as < action_schemas.size(); ++as) {
            const ActionSchema &action_schema = action_schemas[as];
            int64_t size = 1;
            for (int var : action_schema.eff_vars){
                size *= variables[var].get_domain_size();
                if (size > max_leaf_size){
                    break;
                }
            }
            if (size > max_leaf_size){
                continue;
            }

            auto it = leaf_lookup.find(action_schema.eff_vars);
            if (it == leaf_lookup.end()){
                size_t s = potential_leaves.size();
                leaf_lookup[action_schema.eff_vars] = s;
                potential_leaves.emplace_back(action_schema.eff_vars);
                potential_leaves[s].add_leaf_only_schema(as, action_schema);
            } else {
                potential_leaves[it->second].add_leaf_only_schema(as, action_schema);
            }
        }
    }

    log << action_schemas.size() << " action schemes" << endl;
    log << potential_leaves.size() << " potential leaves" << endl;

    vector<vector<size_t>> var_to_p_leaves(task->get_num_variables());
    for (size_t i = 0; i < potential_leaves.size(); ++i) {
        const PotentialLeaf &pleaf = potential_leaves[i];
        for (int var : pleaf.vars) {
            var_to_p_leaves[var].push_back(i);
        }
    }

    // set the number of leaf-only actions
    for (const PotentialLeaf &p_leaf : potential_leaves) {
        vector<size_t> superset_schemes(get_superset_pleaf_ids(p_leaf, potential_leaves, var_to_p_leaves));
        for (size_t superset_as_id : superset_schemes) {
            for (size_t leaf_only_schema : p_leaf.action_schemes) {
                potential_leaves[superset_as_id].add_leaf_only_schema(leaf_only_schema, action_schemas[leaf_only_schema]);
            }
        }
        if (!check_timeout()){
            return;
        }
    }

    if (add_cg_sccs_){
        add_cg_sccs(potential_leaves, var_to_p_leaves);
    }

    vector<vector<unordered_map<size_t, size_t>>> facts_to_mobility;
    vector<vector<size_t>> sum_fact_mobility;
    if (strategy == WMIS_STRATEGY::MFA || min_fact_flexibility > 0) {
        compute_fact_flexibility(facts_to_mobility, sum_fact_mobility);
    }

    multiply_out_action_schemas(potential_leaves, facts_to_mobility, sum_fact_mobility);

    if (!check_timeout()){
        return;
    }

    if (potential_leaves.empty()){
        log << "No potential leaves." << endl;
        return;
    }
}

inline vector<vector<int>> get_sccs(const TaskProxy &task_proxy) {
    size_t num_vars = task_proxy.get_variables().size();
    vector<vector<int> > vars(num_vars);
    const causal_graph::CausalGraph &cg = task_proxy.get_causal_graph();
    for (size_t i = 0; i < num_vars; ++i){
        vars[i] = cg.get_successors(i);
    }
    return sccs::compute_maximal_sccs(vars);
}

void MWISFactoring::add_cg_sccs(vector<PotentialLeaf> &potential_leaves,
                                vector<vector<size_t>> &var_to_p_leaves) {
    vector<vector<int>> sccs = get_sccs(task_proxy);

    if (sccs.size() == 1){
        // we don't want to put all variables into a leaf
        log << "Causal-graph is strongly connected, no potential leaves have been added." << endl;
        return;
    }

    utils::HashSet<vector<int>> leaf_lookup;

    for (const PotentialLeaf &pleaf : potential_leaves){
        leaf_lookup.insert(pleaf.vars);
    }

    size_t added = 0;
    VariablesProxy variables = task_proxy.get_variables();
    for (auto &scc : sccs) {
        int64_t size = 1;
        for (int var : scc){
            size *= variables[var].get_domain_size();
            if (size > max_leaf_size){
                break;
            }
        }
        if (size > max_leaf_size){
            continue;
        }

        sort(scc.begin(), scc.end());

        if (leaf_lookup.find(scc) == leaf_lookup.end()){
            unordered_set<size_t> subset_schemes;
            for (int scc_var : scc) {
                for (size_t index : var_to_p_leaves[scc_var]) {
                    if (index >= potential_leaves.size() - added) {
                        // potential leaves from other SCCs cannot be subset
                        continue;
                    }
                    if (subset_schemes.count(index) > 0){
                        continue;
                    }
                    bool subset_schema = true;
                    for (int var : potential_leaves[index].vars) {
                        if (!binary_search(scc.begin(), scc.end(), var)) {
                            subset_schema = false;
                            break;
                        }
                    }
                    if (subset_schema) {
                        subset_schemes.insert(index);
                    }
                }
            }

            size_t s = potential_leaves.size();
            for (int var : scc) {
                var_to_p_leaves[var].push_back(s);
            }
            potential_leaves.emplace_back(scc);
            ++added;

            for (size_t pleaf : subset_schemes){
                assert(pleaf < potential_leaves.size() - added);
                for (size_t as_id : potential_leaves[pleaf].action_schemes){
                    potential_leaves[s].add_leaf_only_schema(as_id, action_schemas[as_id]);
                }
            }
        }
        if (!check_timeout()){
            return;
        }
    }

    log << "Added " << added << " causal-graph SCC potential leaves." << endl;
}

void MWISFactoring::add_options_to_parser(plugins::Feature &feature) {
    feature.add_option<WMIS_STRATEGY>(
            "strategy",
            "This option determines the property of the factoring that is being "
                  "optimized by the LP, e.g. the number of mobile leaves, or the sum"
                  "of leaf mobility.",
            "MFA"
    );
    feature.add_option<int>(
            "min_mobility",
            "Minimum number of leaf-only actions per leaf factor.",
            "1"
    );
    feature.add_option<double>(
            "min_flexibility",
            "Minimum flexibility (ratio between the number of leaf-only vs. all actions affecting a leaf.",
            "0.2"
    );
    feature.add_option<double>(
            "min_fact_flexibility",
            "Fact flexibility is measured as the mean ratio across all facts in a leaf of the number of"
                  "leaf-only vs. all actions affecting that leaf. This option imposes a minimum on that metric.",
            "0"
    );
    feature.add_option<bool>(
            "add_cg_sccs",
            "If true, every SCC of the causal graph is considered a leaf candidate.",
            "true"
    );
}

static plugins::TypedEnumPlugin<WMIS_STRATEGY> _enum_plugin({
    {"MML", "maximize mobile leaves"},
    {"MMAS", "maximize mobile action schemas"},
    {"MM_OPT", "maximize mobility"},
    {"MFA", "maximize mobile facts"},
    {"MM", "maximize mobility (sum)"},
    {"MCL", "maximize number of mobile conclusive leaves"},
    {"MCM", "maximize conclusive mobility, i.e. number of conclusive actions (sum)"},
});

class MWISFactoringFeature : public plugins::TypedFeature<Factoring, MWISFactoring> {
public:
    MWISFactoringFeature() : TypedFeature("wmis") {
        document_title("Maximum-weight independent set factoring");

        Factoring::add_options_to_feature(*this);
        MWISFactoring::add_options_to_parser(*this);
    }
};

static plugins::FeaturePlugin<MWISFactoringFeature> _plugin;
}

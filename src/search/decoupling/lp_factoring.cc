#include "lp_factoring.h"

#include "../algorithms/sccs.h"

#include "../lp/lp_solver.h"

#include "../plugins/plugin.h"

#include "../utils/timer.h"
#include "../utils/hash.h"

#include "../task_utils/causal_graph.h"
#include "../task_proxy.h"
#include "../tasks/root_task.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>


using namespace std;


namespace decoupling {
void LPFactoring::PotentialLeaf::add_leaf_only_schema(size_t action_schema) {
    if (find(action_schemes.begin(),
             action_schemes.end(),
             action_schema) == action_schemes.end()) {
        action_schemes.push_back(action_schema);
        num_actions += factoring->action_schemas[action_schema].num_actions;
        bool all_in = true;
        for (int pre_var : factoring->action_schemas[action_schema].pre_vars) {
            if (!binary_search(vars.begin(), vars.end(), pre_var)) {
                all_in = false;
                break;
            }
        }
        if (all_in) {
            self_mobile_as.push_back(action_schema);
        }
    }
}

LPFactoring::LPFactoring(const plugins::Options &opts) : Factoring(opts),
                                                         strategy(opts.get<STRATEGY>("strategy")),
                                                         min_mobility(opts.get<int>("min_mobility")),
                                                         min_flexibility(opts.get<double>("min_flexibility")),
                                                         min_fact_flexibility(opts.get<double>("min_fact_flexibility")),
                                                         add_cg_sccs_(opts.get<bool>("add_cg_sccs")) {
    if (log.is_at_least_normal()) {
        log << "Using LP factoring with strategy: ";
        switch (strategy) {
        case STRATEGY::MML:
            log << "maximize number of mobile leaves." << endl;
            break;
        case STRATEGY::MMAS:
            log << "maximize number of mobile action schemas." << endl;
            break;
        case STRATEGY::MM_OPT:
            log << "maximize leaf mobility (exact)." << endl;
            break;
        case STRATEGY::MM_APPROX:
            log << "maximize leaf mobility (approximate)." << endl;
            break;
        case STRATEGY::MFA:
            log << "maximize number of mobile facts." << endl;
            break;
        case STRATEGY::MM:
            log << "maximize leaf mobility (sum)." << endl;
            break;
        case STRATEGY::MCL:
            log << "maximize number of mobile conclusive leaves." << endl;
            break;
        case STRATEGY::MCM:
            log << "maximize conclusive leaf mobility." << endl;
            break;
        default:
            log << "ERROR: unknown LP factoring strategy." << endl;
            exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
        }
    }

    if (strategy != STRATEGY::MFA && min_fact_flexibility > 0.0) {
        log << "Option min_fact_flexibility is only possible in combination with strategy MFA." << endl;
        exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }
}

inline void get_combinations(const vector<size_t> &act_schemas,
                             vector<vector<size_t>> &combinations,
                             vector<size_t> &current,
                             size_t i = 0) {
    if (i == act_schemas.size()) {
        combinations.push_back(current);
        current.pop_back();
        return;
    }
    for (size_t as_num = i; as_num < act_schemas.size(); ++as_num) {
        current.push_back(act_schemas[as_num]);
        get_combinations(act_schemas, combinations, current, as_num + 1);
    }
    if (i != 0) {
        combinations.push_back(current);
        current.pop_back();
    }
}

inline double get_log(size_t num_actions) {
    return log(max(1.0001, (double)num_actions));
}

vector<size_t> LPFactoring::add_center_variables_and_get_ids(named_vector::NamedVector<lp::LPVariable> &variables,
                                                             named_vector::NamedVector<lp::LPConstraint> &constraints,
                                                             const vector<bool> &can_be_leaf_var) const {
    vector<size_t> c_vars_ids(task->get_num_variables(), -1);
    // binary var for each FDR variable; 1 => is in center, 0 => not
    for (int var = 0; var < task->get_num_variables(); ++var) {
        if (can_be_leaf_var[var]) {
            c_vars_ids[var] = variables.size();
            variables.push_back(lp::LPVariable(0.0, 1.0, 0.0, true));
        }
    }

    // set center vars if potential leaf is not a leaf
    for (const auto &pleaf : potential_leaves) {
        // only consider the variables that can possibly be leaf variables
        int num_vars = 0;
        for (int var : pleaf.vars) {
            if (can_be_leaf_var[var]) {
                ++num_vars;
            }
        }
        lp::LPConstraint constraint(-infty, num_vars);
        constraint.insert(pleaf.id, num_vars);
        // if pleaf becomes leaf, then none of its variables can be center vars
        for (int var : pleaf.vars) {
            if (can_be_leaf_var[var]) {
                constraint.insert(c_vars_ids[var], 1.0);
            }
        }
        constraints.push_back(constraint);
    }

    for (int var = 0; var < task->get_num_variables(); ++var) {
        if (can_be_leaf_var[var]) {
            // force LP var to be 1 (var is in center) if no potential leaf that contains var becomes a leaf
            lp::LPConstraint constraint(-infty, -1.0);
            constraint.insert(c_vars_ids[var], -1.0);
            for (auto pleaf_id: var_to_p_leaves[var]) {
                constraint.insert(pleaf_id, -1.0);
            }
            constraints.push_back(constraint);
        }
    }

    return c_vars_ids;
}

void LPFactoring::add_leaf_intersection_constraints(named_vector::NamedVector<lp::LPConstraint> &constraints) const {
    // non-empty intersection between potential leaves
    vector<vector<bool>> pleaf_intersect(potential_leaves.size() - 1);
    for (size_t i = 0; i < pleaf_intersect.size(); ++i) {
        pleaf_intersect[i].resize(potential_leaves.size() - i - 1, false);
    }
    for (int var = 0; var < (int)task->get_num_variables(); ++var) {
        for (int pot_leaf_1 : var_to_p_leaves[var]) {
            for (int pot_leaf_2 : var_to_p_leaves[var]) {
                if (pot_leaf_1 < pot_leaf_2) {
                    pleaf_intersect[pot_leaf_1][pot_leaf_2 - pot_leaf_1 - 1] = true;
                }
            }
        }
    }
    for (size_t p_leaf_1 = 0; p_leaf_1 < potential_leaves.size(); ++p_leaf_1) {
        for (size_t p_leaf_2 = p_leaf_1 + 1; p_leaf_2 < potential_leaves.size(); ++p_leaf_2) {
            // we need p_leaf_1 < p_leaf_2
            if (pleaf_intersect[p_leaf_1][p_leaf_2 - p_leaf_1 - 1]) {
                lp::LPConstraint constraint(0.0, 1.0);
                constraint.insert(p_leaf_1, 1.0);
                constraint.insert(p_leaf_2, 1.0);
                constraints.push_back(constraint);
            }
        }
    }
}

void LPFactoring::add_min_flexibility_constraints(named_vector::NamedVector<lp::LPConstraint> &constraints,
                                                  const vector<vector<size_t>> &mob_as_var_ids) const {
    for (const PotentialLeaf &pleaf : potential_leaves) {
        lp::LPConstraint constraint(-infty, 0.0);
        constraint.insert(pleaf.id, min_flexibility);
        for (size_t as_num = 0; as_num < pleaf.action_schemes.size(); ++as_num) {
            constraint.insert(mob_as_var_ids[pleaf.id][as_num], -pleaf.as_flexibility[as_num]);
        }
        constraints.push_back(constraint);
    }
}

void LPFactoring::add_min_mobility_constraints(named_vector::NamedVector<lp::LPConstraint> &constraints,
                                               const vector<vector<size_t>> &mob_as_var_ids) const {
    for (const PotentialLeaf &pleaf : potential_leaves) {
        lp::LPConstraint constraint(-infty, 0.0);
        constraint.insert(pleaf.id, min_mobility);
        for (size_t as_num = 0; as_num < pleaf.action_schemes.size(); ++as_num) {
            constraint.insert(mob_as_var_ids[pleaf.id][as_num], -action_schemas[pleaf.action_schemes[as_num]].num_actions);
        }
        constraints.push_back(constraint);
    }
}

void LPFactoring::add_potential_leaf_to_action_schema_constraints(named_vector::NamedVector<lp::LPConstraint> &constraints,
                                                                  const vector<vector<size_t>> &mob_as_var_ids,
                                                                  const vector<bool> &can_be_leaf_var,
                                                                  const vector<size_t> &c_vars_ids) const {
    bool skip_self_mobile_leaves = strategy == STRATEGY::MML &&
        min_flexibility == 0.0 &&
        min_mobility <= 1;

    // at least one action schema needs to be mobile in each leaf
    for (const auto &pleaf : potential_leaves) {
        if (skip_self_mobile_leaves && !pleaf.self_mobile_as.empty()) {
            continue;
        }
        lp::LPConstraint constraint(-infty, 0.0);
        constraint.insert(pleaf.id, 1.0);
        for (size_t as_num = 0; as_num < pleaf.action_schemes.size(); ++as_num) {
            constraint.insert(mob_as_var_ids[pleaf.id][as_num], -1.0);
        }
        constraints.push_back(constraint);
    }

    if (!check_timeout()) {
        return;
    }

    // set action schemas mobile if (1) it is part of a leaf, and
    // (2) precondition variables outside the leaf are in the center
    for (const auto &pleaf : potential_leaves) {
        if (skip_self_mobile_leaves && !pleaf.self_mobile_as.empty()) {
            continue;
        }
        for (size_t as_num = 0; as_num < pleaf.action_schemes.size(); ++as_num) {
            vector<int> outside_pre_vars;
            for (int var : action_schemas[pleaf.action_schemes[as_num]].pre_vars) {
                if (can_be_leaf_var[var] &&
                    !binary_search(pleaf.vars.begin(),
                                   pleaf.vars.end(),
                                   var)) {
                    outside_pre_vars.push_back(var);
                }
            }

            // action schema is mobile only if it is a leaf and all pre variables are center
            lp::LPConstraint constraint(-infty, 0.0);
            constraint.insert(pleaf.id, -1.0);
            constraint.insert(mob_as_var_ids[pleaf.id][as_num], outside_pre_vars.size() + 1.0);
            for (int var: outside_pre_vars) {
                constraint.insert(c_vars_ids[var], -1.0);
            }
            constraints.push_back(constraint);
        }
    }

    // an action schema can only be mobile for a leaf if it actually is a leaf
    for (const auto &pleaf : potential_leaves) {
        if (skip_self_mobile_leaves && !pleaf.self_mobile_as.empty()) {
            continue;
        }
        for (size_t as_num = 0; as_num < pleaf.action_schemes.size(); ++as_num) {
            lp::LPConstraint constraint(-infty, 0.0);
            constraint.insert(pleaf.id, -1.0);
            constraint.insert(mob_as_var_ids[pleaf.id][as_num], 1.0);
            constraints.push_back(constraint);
        }
    }
}

bool LPFactoring::is_as_leaf_irrelevant(const ActionSchema &as, const PotentialLeaf &leaf) const {
    assert(!has_as_pre_or_eff_on_leaf(as, leaf));
    vector<bool> is_leaf_pre_var(task->get_num_variables(), false);
    for (auto as_id : leaf.action_schemes) {
        for (int var : action_schemas[as_id].pre_vars) {
            is_leaf_pre_var[var] = true;
        }
    }
    for (int evar : as.eff_vars) {
        if (is_leaf_pre_var[evar]) {
            return false;
        }
    }
    return true;
}

bool LPFactoring::is_as_leaf_conclusive(const ActionSchema &as, const PotentialLeaf &leaf) const {
    // is every var in leaf contained in either the precondition of effect of as
    for (int lvar : leaf.vars) {
        bool covered = false;
        for (int pvar: as.pre_vars) {
            if (lvar == pvar) {
                covered = true;
            } else if (pvar > lvar) {
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
            if (!covered) {
                return false;
            }
        }
    }
    return true;
}

bool LPFactoring::has_as_pre_or_eff_on_leaf(const ActionSchema &as, const PotentialLeaf &leaf) const {
    for (int lvar : leaf.vars) {
        for (int pvar: as.pre_vars) {
            if (lvar == pvar) {
                return true;
            } else if (pvar > lvar) {
                break;
            }
        }
        for (int evar: as.eff_vars) {
            if (lvar == evar) {
                return true;
            } else if (evar > lvar) {
                break;
            }
        }
    }
    return false;
}

void LPFactoring::construct_lp_conclusive(named_vector::NamedVector<lp::LPVariable> &variables,
                                          named_vector::NamedVector<lp::LPConstraint> &constraints) {
    assert(variables.size() == 0);
    assert(constraints.size() == 0);

    compute_action_schemas();

    if (action_schemas.empty()) {
        // mostly for trivially unsolvable task from translator?
        log << "ERROR: No action schemas." << endl;
        return;
    }

    compute_potential_leaves();

    if (static_cast<int>(potential_leaves.size()) < min_number_leaves) {
        log << "Only " << potential_leaves.size() <<
            " potential leaves left, but minimum number of leaves is " <<
            min_number_leaves << "." << endl;
        return;
    }

    if (!check_timeout()) {
        return;
    }

    vector<bool> can_be_leaf_var(task->get_num_variables(), false);
    // binary var for each potential leaf; 1 => becomes leaf, 0 => not
    for (const auto &pleaf : potential_leaves) {
        for (int var : pleaf.vars) {
            can_be_leaf_var[var] = true;
        }
        double obj = 0.0;
        if (strategy == STRATEGY::MCL) {
            // for tie-breaking, give objective value of 1 to every mobile leaf
            obj = 1.0;
        }
        variables.push_back(lp::LPVariable(0.0, 1.0, obj, true));
    }

    // an LP variable for every potential leaf that is 1 iff the leaf is conclusive
    vector<int> concl_pleaf_lp_vars;
    if (strategy == STRATEGY::MCL) {
        assert(variables.size() == static_cast<int>(potential_leaves.size()));
        for (size_t i = 0; i < potential_leaves.size(); ++i) {
            concl_pleaf_lp_vars.push_back(variables.size());
            variables.push_back(lp::LPVariable(0.0, 1.0, 100.0, true));
        }
        assert(concl_pleaf_lp_vars.size() == potential_leaves.size());
        assert(variables.size() == 2 * static_cast<int>(potential_leaves.size()));
    }

    // binary var for each action schema of each potential leaf; 1 => is mobile, 0 => not
    vector<vector<size_t>> mob_as_var_ids(potential_leaves.size());
    // binary var for each action schema of each potential leaf; 1 => is conclusive for potential leaf, 0 => not
    vector<vector<size_t>> concl_as_var_ids(potential_leaves.size());

    // for every potential leaf we store the list of action schemas relevant for the constraints;
    // for MCL, this is the schemas that can be neither irrelevant nor conclusive for the leaf,
    // for MCM, it's the schemas that become irrelevant or conclusive as global actions
    vector<vector<size_t>> relevant_as_by_pleaf(potential_leaves.size());

    // list of LP variable IDs that represent mob_as_var_ids for every action schemas, i.e.,
    // the corresponding LP variable IDs across potential leaves
    vector<vector<size_t>> as_to_mob_as_var_ids(action_schemas.size());

    for (const PotentialLeaf &pleaf : potential_leaves) {
        for (auto as_id : pleaf.action_schemes) {
            as_to_mob_as_var_ids[as_id].push_back(variables.size());
            mob_as_var_ids[pleaf.id].push_back(variables.size());
            variables.push_back(lp::LPVariable(0.0, 1.0, 0.0, true));
        }
        // for all action schemas as_id that are not part of pleaf, add a variable that represents if as_id is conclusive for pleaf
        for (size_t as_id = 0; as_id < action_schemas.size(); ++as_id) {
            if (find(pleaf.action_schemes.begin(), pleaf.action_schemes.end(), as_id) != pleaf.action_schemes.end()) {
                // action schema is in pleaf => cannot be conclusive for pleaf
                continue;
            }
            bool conclusive_or_irrelevant = false;
            if (has_as_pre_or_eff_on_leaf(action_schemas[as_id], pleaf)) {
                if (is_as_leaf_conclusive(action_schemas[as_id], pleaf)) {
                    conclusive_or_irrelevant = true;
                }
            } else if (is_as_leaf_irrelevant(action_schemas[as_id], pleaf)) {
                conclusive_or_irrelevant = true;
            }
            if (strategy == STRATEGY::MCM) {
                if (!conclusive_or_irrelevant) {
                    // the action schema cannot be conclusive for this potential leaf, no need to add variables for the objective value
                    continue;
                }
            } else {
                assert(strategy == STRATEGY::MCL);
                if (conclusive_or_irrelevant) {
                    // the action schema will be conclusive or irrelevant for the leaf, no need to add a constraint
                    continue;
                }
            }

            relevant_as_by_pleaf[pleaf.id].push_back(as_id);

            if (strategy == STRATEGY::MCM) {
                concl_as_var_ids[pleaf.id].push_back(variables.size());
                variables.push_back(lp::LPVariable(0.0,
                                                   1.0,
                                                   action_schemas[as_id].num_actions,
                                                   true));
            }
        }
    }

    if (min_flexibility > 0.0) {
        add_min_flexibility_constraints(constraints, mob_as_var_ids);
    }

    if (min_mobility > 1) {
        add_min_mobility_constraints(constraints, mob_as_var_ids);
    }

    vector<size_t> c_vars_ids = add_center_variables_and_get_ids(variables, constraints, can_be_leaf_var);

    add_leaf_intersection_constraints(constraints);

    if (!check_timeout()) {
        return;
    }

    add_potential_leaf_to_action_schema_constraints(constraints, mob_as_var_ids, can_be_leaf_var, c_vars_ids);

    // add constraints such that conclusiveness of action schemas for pleaves is set properly
    for (const auto &pleaf : potential_leaves) {
        if (strategy == STRATEGY::MCL) {
            lp::LPConstraint constraint1(-infty, 0.0);
            constraint1.insert(pleaf.id, -1.0);
            constraint1.insert(concl_pleaf_lp_vars[pleaf.id], 1.0);
            // pleaf can only be conclusive if it becomes a (mobile) leaf
            constraints.push_back(constraint1);
        }
        for (size_t i = 0; i < relevant_as_by_pleaf[pleaf.id].size(); ++i) {
            size_t as_id = relevant_as_by_pleaf[pleaf.id][i];
            if (strategy == STRATEGY::MCM) {
                lp::LPConstraint constraint1(-infty, 0.0);
                constraint1.insert(pleaf.id, -1.0);
                constraint1.insert(concl_as_var_ids[pleaf.id][i], 1.0);
                // action schema can only be conclusive for pleaf if pleaf becomes a leaf
                constraints.push_back(constraint1);

                // action schema can be conclusive for pleaf only if it is a global action schema, i.e., iff it is not mobile for any potential leaf
                for (auto as_mob_var_id : as_to_mob_as_var_ids[as_id]) {
                    lp::LPConstraint constraint2(-infty, 1.0);
                    constraint2.insert(concl_as_var_ids[pleaf.id][i], 1.0);
                    constraint2.insert(as_mob_var_id, 1.0);
                    constraints.push_back(constraint2);
                }
            } else {
                assert(strategy == STRATEGY::MCL);

                // action schema cannot be conclusive or irrelevant for pleaf so needs to be leaf-only for some potential leaf
                lp::LPConstraint constraint2(-infty, 0.0);
                constraint2.insert(concl_pleaf_lp_vars[pleaf.id], 1.0);
                for (auto as_mob_var_id : as_to_mob_as_var_ids[as_id]) {
                    constraint2.insert(as_mob_var_id, -1.0);
                }
                constraints.push_back(constraint2);
            }
        }
    }
}

void LPFactoring::construct_lp_all(named_vector::NamedVector<lp::LPVariable> &variables,
                                   named_vector::NamedVector<lp::LPConstraint> &constraints) {
    assert(variables.size() == 0);
    assert(constraints.size() == 0);

    VariablesProxy vars_proxy = task_proxy.get_variables();

    compute_action_schemas();

    vector<vector<unordered_map<size_t, size_t>>> facts_to_mobility;
    vector<vector<size_t>> sum_fact_mobility;
    if (strategy == STRATEGY::MFA) {
        // need to store information for effect *facts*, which gets otherwise lost
        facts_to_mobility.resize(task->get_num_variables());
        for (int var = 0; var < task->get_num_variables(); ++var) {
            facts_to_mobility[var].resize(vars_proxy[var].get_domain_size());
        }

        utils::HashMap<vector<int>, utils::HashMap<vector<int>, size_t>> scheme_lookup;
        for (size_t as_id = 0; as_id < action_schemas.size(); ++as_id) {
            const ActionSchema &as = action_schemas[as_id];
            scheme_lookup[as.pre_vars][as.eff_vars] = as_id;
        }

        for (OperatorProxy op : task_proxy.get_operators()) {
            vector<int> pre_vars;
            for (FactProxy pre : op.get_preconditions()) {
                pre_vars.push_back(pre.get_variable().get_id());
            }
            sort(pre_vars.begin(), pre_vars.end());

            vector<int> eff_vars;
            for (EffectProxy eff : op.get_effects()) {
                eff_vars.push_back(eff.get_fact().get_variable().get_id());
            }
            sort(eff_vars.begin(), eff_vars.end());

            assert(scheme_lookup.count(pre_vars) > 0 && scheme_lookup[pre_vars].count(eff_vars) > 0);
            size_t as = scheme_lookup[pre_vars][eff_vars];
            for (EffectProxy eff : op.get_effects()) {
                facts_to_mobility[eff.get_fact().get_variable().get_id()][eff.get_fact().get_value()][as]++;
            }
        }

        sum_fact_mobility.resize(task->get_num_variables());
        for (int var = 0; var < task->get_num_variables(); ++var) {
            sum_fact_mobility[var].resize(vars_proxy[var].get_domain_size());
            for (int val = 0; val < vars_proxy[var].get_domain_size(); ++val) {
                for (const auto & [as, num] : facts_to_mobility[var][val]) {
                    sum_fact_mobility[var][val] += num;
                }
            }
        }
    }

    if (action_schemas.empty()) {
        // mostly for trivially unsolvable task from translator?
        log << "ERROR: No action schemas." << endl;
        return;
    }

    compute_potential_leaves();

    if (!check_timeout()) {
        return;
    }

    if (static_cast<int>(potential_leaves.size()) < min_number_leaves) {
        log << "Only " << potential_leaves.size() <<
            " potential leaves left, but minimum number of leaves is " <<
            min_number_leaves << "." << endl;
        return;
    }

    // LP variables
    vector<double> min_leaf_mobility;
    vector<size_t> non_self_mobile_leaf_actions;
    if (strategy == STRATEGY::MM_APPROX) {
        min_leaf_mobility.resize(potential_leaves.size(), 0.0);
        non_self_mobile_leaf_actions.resize(potential_leaves.size(), 0);
    }
    vector<bool> can_be_leaf_var(task->get_num_variables(), false);
    // binary var for each potential leaf; 1 => becomes leaf, 0 => not
    for (const auto &pleaf : potential_leaves) {
        double obj = 0.0;
        if (strategy == STRATEGY::MM_APPROX) {
            // here obj is minimum leaf mobility
            non_self_mobile_leaf_actions[pleaf.id] = pleaf.num_actions;
            if (!pleaf.self_mobile_as.empty()) {
                size_t mob = 0;
                for (size_t act_schema : pleaf.self_mobile_as) {
                    mob += action_schemas[act_schema].num_actions;
                    non_self_mobile_leaf_actions[pleaf.id] -= action_schemas[act_schema].num_actions;
                }
                obj = get_log(mob);
                min_leaf_mobility[pleaf.id] = obj;
            }
        } else if (strategy == STRATEGY::MML) {
            // here, obj is 1 for every potential leaf
            obj = 1.0;
        }
        for (int var : pleaf.vars) {
            can_be_leaf_var[var] = true;
        }
        variables.push_back(lp::LPVariable(0.0, 1.0, obj, true));
    }
    size_t act_schema_offset = variables.size();
    if (strategy == STRATEGY::MFA) {
        // binary variable for each action schema: 1 -> mobile, 0 -> not
        for (size_t as = 0; as < action_schemas.size(); ++as) {
            variables.push_back(lp::LPVariable(0.0, 1.0, 0.0, true));
        }
    }
    // binary var for each action schema of each potential leaf; 1 => is mobile, 0 => not
    vector<vector<size_t>> mob_as_var_ids(potential_leaves.size());
    vector<vector<size_t>> as_pleaf_var_ids;
    if (strategy == STRATEGY::MFA) {
        // ids to the LP variables where an action schema is part of a potential leaf
        as_pleaf_var_ids.resize(action_schemas.size());
    }
    for (const PotentialLeaf &pleaf : potential_leaves) {
        if (strategy == STRATEGY::MML &&
            !pleaf.self_mobile_as.empty() &&
            min_flexibility == 0.0 &&
            min_mobility <= 1) {
            // if the leaf is self-mobile, there is no need to force
            // one of its action schemas to be mobile
            continue;
        }
        for (size_t act_schema : pleaf.action_schemes) {
            mob_as_var_ids[pleaf.id].push_back(variables.size());
            switch (strategy) {
            case STRATEGY::MMAS:
                variables.push_back(lp::LPVariable(0.0, 1.0, 1.0, true));
                break;
            case STRATEGY::MFA:
                as_pleaf_var_ids[act_schema].push_back(variables.size());
            // fall through
            case STRATEGY::MML:     // fall through
            case STRATEGY::MM_OPT:
                variables.push_back(lp::LPVariable(0.0, 1.0, 0.0, true));
                break;
            case STRATEGY::MM:
                variables.push_back(lp::LPVariable(0.0, 1.0, action_schemas[act_schema].num_actions, true));
                break;
            case STRATEGY::MM_APPROX:
            {
                double obj_v = 0.0;
                if (find(pleaf.self_mobile_as.begin(), pleaf.self_mobile_as.end(), act_schema) == pleaf.self_mobile_as.end()) {
                    // otherwise, mobility is accounted for by the objective value of pleaf.id
                    obj_v = action_schemas[act_schema].num_actions *
                        (get_log(pleaf.num_actions) - min_leaf_mobility[pleaf.id]) /
                        non_self_mobile_leaf_actions[pleaf.id];
                }
                variables.push_back(lp::LPVariable(0.0, 1.0, obj_v, true));
            }
            break;
            default:
                log << "unknown strategy" << endl;
                exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
            }
        }
    }
    vector<vector<size_t>> as_combs_ids;
    if (strategy == STRATEGY::MM_OPT) {
        as_combs_ids.resize(potential_leaves.size());
        // binary var for each combination of mobile action schemas of each potential leaf; 1 => combination is mobile, 0 => not
        for (const PotentialLeaf &pleaf : potential_leaves) {
            vector<vector<size_t>> combinations;
            vector<size_t> tmp;
            get_combinations(pleaf.action_schemes, combinations, tmp);
            for (const auto &act_schemas : combinations) {
                as_combs_ids[pleaf.id].push_back(variables.size());
                size_t num_actions = 0;
                for (size_t as : act_schemas) {
                    num_actions += action_schemas[as].num_actions;
                }
                variables.push_back(lp::LPVariable(0.0, 1.0, get_log(num_actions), true));
            }
        }
    }

    vector<size_t> c_vars_ids = add_center_variables_and_get_ids(variables, constraints, can_be_leaf_var);

    vector<size_t> mob_cvars_ids;
    if (strategy == STRATEGY::MFA) {
        // float var for each FDR fact; 1 => fully mobile, 0 => not mobile
        mob_cvars_ids.resize(task->get_num_variables());
        for (int var = 0; var < task->get_num_variables(); ++var) {
            if (can_be_leaf_var[var]) {
                mob_cvars_ids[var] = variables.size();
                for (int val = 0; val < vars_proxy[var].get_domain_size(); ++val) {
                    variables.push_back(lp::LPVariable(0.0, 1.0, 1.0));
                }
            }
        }
    }

    // LP constraints
    if (min_flexibility > 0.0) {
        add_min_flexibility_constraints(constraints, mob_as_var_ids);
    }

    if (min_mobility > 1) {
        add_min_mobility_constraints(constraints, mob_as_var_ids);
    }

    if (!check_timeout()) {
        return;
    }

    add_leaf_intersection_constraints(constraints);

    if (!check_timeout()) {
        return;
    }

    add_potential_leaf_to_action_schema_constraints(constraints, mob_as_var_ids, can_be_leaf_var, c_vars_ids);

    if (!check_timeout()) {
        return;
    }

    if (strategy == STRATEGY::MM_OPT) {
        // for each potential leaf, exactly one combination of its action schemas is mobile
        for (const PotentialLeaf &pleaf : potential_leaves) {
            size_t i = 0;
            vector<vector<size_t>> combinations;
            vector<size_t> tmp;
            get_combinations(pleaf.action_schemes, combinations, tmp);
            for (const auto &c : combinations) {
                lp::LPConstraint constraint(-infty, 0.0);
                constraint.insert(pleaf.id, -1.0);
                constraint.insert(as_combs_ids[pleaf.id][i], 1.0 + c.size());
                for (size_t as : c) {
                    // TODO build a map for this
                    size_t id = find(pleaf.action_schemes.begin(), pleaf.action_schemes.end(), as) - pleaf.action_schemes.begin();
                    constraint.insert(mob_as_var_ids[pleaf.id][id], -1.0);
                }
                constraints.push_back(constraint);

                if (c.size() < pleaf.action_schemes.size()) {
                    lp::LPConstraint constraint1(-infty, pleaf.action_schemes.size() - c.size());
                    constraint1.insert(as_combs_ids[pleaf.id][i], pleaf.action_schemes.size() - c.size());
                    for (size_t as_num = 0; as_num < pleaf.action_schemes.size(); ++as_num) {
                        // TODO avoid searching
                        if (find(c.begin(), c.end(), pleaf.action_schemes[as_num]) == c.end()) {
                            constraint1.insert(mob_as_var_ids[pleaf.id][as_num], 1.0);
                        }
                    }
                    constraints.push_back(constraint1);
                }
                ++i;
            }
        }
    }

    if (!check_timeout()) {
        return;
    }

    if (strategy == STRATEGY::MFA) {
        // an action schema is mobile if it is mobile for a potential leaf
        for (size_t as_id = 0; as_id < action_schemas.size(); ++as_id) {
            lp::LPConstraint constraint(-infty, 0.0);
            constraint.insert(act_schema_offset + as_id, 1.0);
            for (size_t as_mob_id : as_pleaf_var_ids[as_id]) {
                constraint.insert(as_mob_id, -1.0);
            }
            constraints.push_back(constraint);
        }

        // fact mobility is upper bounded by sum of mobility ratios over affecting action schemas
        for (int var = 0; var < task->get_num_variables(); ++var) {
            if (can_be_leaf_var[var]) {
                for (int val = 0; val < vars_proxy[var].get_domain_size(); ++val) {
                    lp::LPConstraint constraint(-infty, 0.0);
                    constraint.insert(mob_cvars_ids[var] + val, 1.0);
                    for (const auto & [as, num] : facts_to_mobility[var][val]) {
                        constraint.insert(act_schema_offset + as, -(double)num / sum_fact_mobility[var][val]);
                    }
                    constraints.push_back(constraint);
                }
            }
        }

        if (min_fact_flexibility > 0.0) {
            // average flexibility of all facts in the leaf must be >= min_fact_flexibility
            for (const PotentialLeaf &pleaf : potential_leaves) {
                lp::LPConstraint constraint(-infty, 0.0);
                int num_facts = 0;
                for (int var : pleaf.vars) {
                    num_facts += vars_proxy[var].get_domain_size();
                    for (int val = 0; val < vars_proxy[var].get_domain_size(); ++val) {
                        constraint.insert(mob_cvars_ids[var] + val, -1.0);
                    }
                }
                constraint.insert(pleaf.id, num_facts * min_fact_flexibility);
                constraints.push_back(constraint);
            }
        }
    }
}

void LPFactoring::compute_factoring_() {
    lp::LPSolver solver(lp::LPSolverType::CPLEX);
    infty = solver.get_infinity();

    named_vector::NamedVector<lp::LPVariable> variables;
    named_vector::NamedVector<lp::LPConstraint> constraints;

    if (strategy == STRATEGY::MCL || strategy == STRATEGY::MCM) {
        construct_lp_conclusive(variables, constraints);
    } else {
        construct_lp_all(variables, constraints);
    }

    if (!check_timeout()) {
        return;
    }

    if (variables.size() == 0) {
        log << "WARNING: no LP variables created, stopping." << endl;
        return;
    }

    // save memory
    vector<ActionSchema>().swap(action_schemas);

    if (min_number_leaves >= 2) {
        // adding the constraint typically increases runtime of CPLEX,
        // therefore, we only add it if we need strictly more than 2 leaves,
        // since 2 leaves are guaranteed by the elimination of potential leaves
        // in set_potential_leaves()
        lp::LPConstraint constraint(min_number_leaves, min(static_cast<int>(potential_leaves.size()), task->get_num_variables()));
        for (size_t p_leaf = 0; p_leaf < potential_leaves.size(); ++p_leaf) {
            constraint.insert(p_leaf, 1.0);
        }
        constraints.push_back(constraint);
    }

    vector<double> solution;
    if (constraints.size() == 0) {
        // this should only happen when all potential leaves can become leaves
        solution.resize(potential_leaves.size(), 1);
        log << "No constraints, all candidates become leaf factors." << endl;
    } else {
        solver.load_problem(lp::LinearProgram(lp::LPObjectiveSense::MAXIMIZE,
                                              move(variables),
                                              move(constraints),
                                              infty));

        solver.set_time_limit(factoring_timer.get_remaining_time());

        solver.solve();

        if (log.is_at_least_verbose()) {
            solver.print_statistics();
            solver.print_failure_analysis();
            if (solver.has_solution()) {
                log << "LP objective value: " << solver.get_objective_value() << endl;
            }
        }

        if (solver.has_solution()) {
            solution = solver.extract_solution();
        }

        if (solution.empty()) {
            log << "WARNING: no solution found by LP solver." << endl;
            return;
        }
    }

    bool leaves_overlap = false;
    vector<bool> is_leaf_var(task->get_num_variables(), false);
    double epsilon = 0.01;
    for (size_t p_leaf = 0; p_leaf < potential_leaves.size(); ++p_leaf) {
        if (static_cast<int>(ceil(solution[p_leaf] - epsilon)) == 1) {
            set<int> leaf;
            for (int var : potential_leaves[p_leaf].vars) {
                if (is_leaf_var[var]) {
                    leaves_overlap = true;
                    // simple sanity check
                }
                is_leaf_var[var] = true;
                leaf.insert(var);
            }
            leaves.emplace_back(leaf.begin(), leaf.end());
        }
    }

    if (leaves_overlap) {
        // TODO: investigate this!
        // this tends to happen when CPLEX is running out of memory on the cluster,
        // when run with an external memory limit. could try to additionally set
        // the memory limit for CPLEX explicitly
        log << "ERROR: solution returned by LP solver is not a valid factoring: leaves do overlap." << endl;
        exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    }
}

void LPFactoring::save_memory() {
    Factoring::save_memory();
    vector<set<int>>().swap(var_to_p_leaves);
    vector<PotentialLeaf>().swap(potential_leaves);
}

void LPFactoring::compute_leaf_flexibility() {
    compute_var_to_ops_map();
    for (auto &pl : potential_leaves) {
        if (!pl.as_flexibility.empty()) {
            continue;
        }
        size_t num_ops = 0;
        vector<bool> handled_op(task->get_num_operators(), false);
        for (int var : pl.vars) {
            for (size_t op_id : var_to_affecting_op[var]) {
                if (!handled_op[op_id]) {
                    handled_op[op_id] = true;
                    ++num_ops;
                }
            }
        }
        pl.as_flexibility = vector<double>(pl.action_schemes.size(), 0.0);
        pl.max_flexibility = 0.0;
        for (size_t as_num = 0; as_num < pl.action_schemes.size(); ++as_num) {
            double as_flex = (double)action_schemas[pl.action_schemes[as_num]].num_actions / num_ops;
            pl.as_flexibility[as_num] = as_flex;
            pl.max_flexibility += as_flex;
        }
    }
}

void LPFactoring::filter_potential_leaves() {
    // filter leaves that cannot possibly have more than the minimum flexibility or mobility
    set<size_t> erase_leaves;
    if (min_flexibility > 0.0) {
        compute_leaf_flexibility();
        for (const auto &pl : potential_leaves) {
            if (pl.max_flexibility < min_flexibility) {
                erase_leaves.insert(pl.id);
            }
        }
    }
    if (min_mobility > 1) {
        for (const auto &pl : potential_leaves) {
            if (pl.num_actions < min_mobility) {
                erase_leaves.insert(pl.id);
            }
        }
    }
    log << "Removed " << erase_leaves.size() << " potential leaves because of minimal flexibility/mobility." << endl;

    for (auto it = erase_leaves.rbegin(); it != erase_leaves.rend(); ++it) {
        assert(potential_leaves[*it].id == *it);
        potential_leaves.erase(potential_leaves.begin() + *it);
    }

    if (!erase_leaves.empty()) {
        size_t id = 0;
        for (auto &pl : potential_leaves) {
            pl.id = id++;
        }
        recompute_var_to_p_leaves();
    }
}

void LPFactoring::recompute_var_to_p_leaves() {
    if (static_cast<int>(var_to_p_leaves.size()) != task->get_num_variables()) {
        var_to_p_leaves = vector<set<int>>(task->get_num_variables(), set<int>());
    }
    for (int var = 0; var < task->get_num_variables(); ++var) {
        var_to_p_leaves[var].clear();
    }
    for (const PotentialLeaf &pleaf : potential_leaves) {
        for (int var : pleaf.vars) {
            var_to_p_leaves[var].insert(pleaf.id);
        }
    }
}

void LPFactoring::compute_potential_leaves() {
    assert(!action_schemas.empty());
    assert(potential_leaves.empty());

    {
        VariablesProxy variables = task_proxy.get_variables();
        utils::HashMap<vector<int>, size_t> leaf_lookup;
        for (size_t as = 0; as < action_schemas.size(); ++as) {
            const ActionSchema &action_schema = action_schemas[as];
            int64_t size = 1;
            for (int var : action_schema.eff_vars) {
                size *= variables[var].get_domain_size();
                if (size > max_leaf_size) {
                    break;
                }
            }
            if (size > max_leaf_size) {
                continue;
            }

            auto it = leaf_lookup.find(action_schema.eff_vars);
            if (it == leaf_lookup.end()) {
                size_t s = potential_leaves.size();
                leaf_lookup[action_schema.eff_vars] = s;
                potential_leaves.emplace_back(this, s, action_schema.eff_vars);
                potential_leaves[s].add_leaf_only_schema(as);
            } else {
                potential_leaves[it->second].add_leaf_only_schema(as);
            }
        }
    }

    log << action_schemas.size() << " action schemes" << endl;
    log << potential_leaves.size() << " potential leaves" << endl;

    recompute_var_to_p_leaves();

    // set the number of leaf-only actions
    for (const PotentialLeaf &pleaf : potential_leaves) {
        if (pleaf.vars.size() == 1) {
            // no need to check for supersets here
            for (size_t index : var_to_p_leaves[pleaf.vars[0]]) {
                for (size_t leaf_only_schema : pleaf.action_schemes) {
                    potential_leaves[index].add_leaf_only_schema(leaf_only_schema);
                }
            }
        } else {
            // since we consider all supersets of p_leaf, we can simply take p_leaf.vars[0],
            // since all supersets must also include that variable
            for (size_t check_pleaf_id : var_to_p_leaves[pleaf.vars[0]]) {
                bool superset_schema = true;
                for (int var : pleaf.vars) {
                    if (!binary_search(potential_leaves[check_pleaf_id].vars.begin(), potential_leaves[check_pleaf_id].vars.end(), var)) {
                        superset_schema = false;
                        break;
                    }
                }
                if (superset_schema) {
                    for (size_t leaf_only_schema : pleaf.action_schemes) {
                        potential_leaves[check_pleaf_id].add_leaf_only_schema(leaf_only_schema);
                    }
                }
            }
        }
        if (!check_timeout()) {
            return;
        }
    }

    if (add_cg_sccs_) {
        add_cg_sccs();
    }

    filter_potential_leaves();

    if (!check_timeout()) {
        return;
    }

    if (potential_leaves.empty()) {
        log << "No potential leaves." << endl;
        return;
    }
}

inline vector<vector<int>> get_sccs(const TaskProxy &task_proxy) {
    size_t num_vars = task_proxy.get_variables().size();
    vector<vector<int>> vars(num_vars);
    const causal_graph::CausalGraph &cg = task_proxy.get_causal_graph();
    for (size_t i = 0; i < num_vars; ++i) {
        vars[i] = cg.get_successors(i);
    }
    return sccs::compute_maximal_sccs(vars);
}

void LPFactoring::add_cg_sccs() {
    vector<vector<int>> sccs = get_sccs(task_proxy);

    if (sccs.size() == 1) {
        // we don't want to put all variables into a leaf
        log << "Causal-graph is strongly connected, no potential leaves have been added." << endl;
        return;
    }

    utils::HashSet<vector<int>> leaf_lookup;

    for (const PotentialLeaf &pleaf : potential_leaves) {
        leaf_lookup.insert(pleaf.vars);
    }

    size_t added = 0;
    VariablesProxy variables = task_proxy.get_variables();
    for (auto &scc : sccs) {
        int64_t size = 1;
        for (int var : scc) {
            size *= variables[var].get_domain_size();
            if (size > max_leaf_size) {
                break;
            }
        }
        if (size > max_leaf_size) {
            continue;
        }

        sort(scc.begin(), scc.end());

        if (leaf_lookup.count(scc) == 0) {
            unordered_set<size_t> subset_schemes;
            for (int scc_var : scc) {
                for (size_t index : var_to_p_leaves[scc_var]) {
                    if (index >= potential_leaves.size() - added) {
                        // potential leaves from other SCCs cannot be subset
                        continue;
                    }
                    if (subset_schemes.count(index) > 0) {
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
            if (subset_schemes.empty()) {
                // this can happen when an SCC is not affected by any action
                continue;
            }

            size_t s = potential_leaves.size();
            for (int var : scc) {
                var_to_p_leaves[var].insert(s);
            }
            potential_leaves.emplace_back(this, s, scc);
            ++added;

            for (size_t pleaf : subset_schemes) {
                assert(pleaf < potential_leaves.size() - added);
                for (size_t as_id : potential_leaves[pleaf].action_schemes) {
                    potential_leaves[s].add_leaf_only_schema(as_id);
                }
            }
        }
        if (!check_timeout()) {
            return;
        }
    }

    log << "Added " << added << " causal-graph SCC potential leaves." << endl;
}

void LPFactoring::add_options_to_parser(plugins::Feature &feature) {
    feature.add_option<STRATEGY>(
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

static plugins::TypedEnumPlugin<STRATEGY> _enum_plugin({
        {"MML", "maximize mobile leaves"},
        {"MMAS", "maximize mobile action schemas"},
        {"MM_OPT", "maximize mobility"},
        {"MM_APPROX", "maximize mobility (approximation)"},
        {"MFA", "maximize mobile facts"},
        {"MM", "maximize mobility (sum)"},
        {"MCL", "maximize number of mobile conclusive leaves"},
        {"MCM", "maximize conclusive mobility, i.e. number of conclusive actions (sum)"},
    });

class LPFactoringFeature : public plugins::TypedFeature<Factoring, LPFactoring> {
public:
    LPFactoringFeature() : TypedFeature("lp") {
        document_title("LP factoring");

        Factoring::add_options_to_feature(*this);
        LPFactoring::add_options_to_parser(*this);
    }
};

static plugins::FeaturePlugin<LPFactoringFeature> _plugin;
}

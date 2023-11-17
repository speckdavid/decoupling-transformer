#include "decoupled_root_task.h"

#include "../axioms.h"
#include "../decoupling/factoring.h"
#include "../operator_id.h"
#include "../plugins/plugin.h"
#include "../task_utils/task_dump.h"
#include "../task_utils/task_properties.h"
#include "../utils/rng_options.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>

using namespace std;

namespace tasks {
DecoupledRootTask::DecoupledRootTask(const plugins::Options &options)
    : RootTask(),
      original_root_task(dynamic_pointer_cast<RootTask>(tasks::g_root_task)),
      rng(utils::parse_rng_from_options(options)),
      factoring(options.get<shared_ptr<decoupling::Factoring>>("factoring")),
      same_leaf_preconditons_single_variable(options.get<bool>("same_leaf_preconditons_single_variable")),
      implicit_effects(options.get<bool>("implicit_effects")) {
    if (implicit_effects) {
        cerr << "not implemented yet!" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_UNSUPPORTED);
    }

    TaskProxy original_task_proxy(*original_root_task);
    task_properties::verify_no_axioms(original_task_proxy);
    task_properties::verify_no_conditional_effects(original_task_proxy);

    factoring->compute_factoring();

    utils::Timer transformation_timer;

    utils::g_log << "Starting decoupling transformaton!" << endl;
    utils::g_log << "Creating new variables..." << endl;
    create_variables();
    utils::g_log << "Creating new mutexes..." << endl;
    create_mutexes();
    utils::g_log << "Creating new operators..." << endl;
    create_operators();
    utils::g_log << "Creating new initial state..." << endl;
    create_initial_state();
    utils::g_log << "Creating new goals..." << endl;
    create_goal();
    utils::g_log << "Creating new axioms..." << endl;
    create_axioms();

    TaskProxy task_proxy(*this);

    // This is also done in the root task which is honestly quite hacky!
    AxiomEvaluator &axiom_evaluator = g_axiom_evaluators[task_proxy];
    axiom_evaluator.evaluate(initial_state_values);
    assert(are_initial_states_consistent());

    utils::g_log << "Time for decoupled transformation: " << transformation_timer << endl;
    print_statistics();

    release_memory();

    if (options.get<bool>("dump_task")) {
        dump();
    }

    if (options.get<bool>("write_sas_file")) {
        write_sas_file("dec_output.sas");
        utils::exit_with(utils::ExitCode::SEARCH_UNSOLVED_INCOMPLETE);
    }
}

int DecoupledRootTask::get_original_operator_id(int op_id) const {
    if (original_op_id_to_global_op_id.count(op_id) == 0) {
        assert(!factoring->is_global_operator(op_id));
        return -1;
    }

    return original_op_id_to_global_op_id.at(op_id);
}

void DecoupledRootTask::reconstruct_plan_if_necessary(vector<OperatorID> &path,
                                                      vector<State> &states) const {
    // remap operator IDs to original operator IDs
    vector<OperatorID> mapped_path;
    for (auto op_id : path) {
        mapped_path.emplace_back(global_op_id_to_original_op_id.at(op_id.get_index()));
    }
    path = mapped_path;
    factoring->insert_leaf_paths(path, states, original_root_task);
}

void DecoupledRootTask::print_statistics() const {
    int num_primary_vars = count_if(variables.begin(), variables.end(), [](const auto &var)
                                    {return var.axiom_layer == -1;});
    int num_secondary_vars = variables.size() - num_primary_vars;

    // utils::g_log << "Original task size: " << original_root_task->get_encoding_size(false) << endl;
    utils::g_log << "Original number of variables: " << original_root_task->variables.size() << endl;
    utils::g_log << "Original number of primary variables: " << original_root_task->variables.size() << endl;
    utils::g_log << "Original number of secondary variables: " << 0 << endl;
    utils::g_log << "Original number of operators: " << get_num_operators() << endl;
    utils::g_log << "Original number of axioms: " << get_num_axioms() << endl;

    // utils::g_log << "Task size: " << get_encoding_size(false) << endl;
    utils::g_log << "Number of variables: " << variables.size() << endl;
    utils::g_log << "Number of primary variables: " << num_primary_vars << endl;
    utils::g_log << "Number of secondary variables: " << num_secondary_vars << endl;
    utils::g_log << "Number of operators: " << get_num_operators() << endl;
    utils::g_log << "Number of axioms: " << get_num_axioms() << endl;
}

void DecoupledRootTask::write_sas_file(const std::string file_name) const {
    utils::Timer write_sas_file_timer;
    utils::g_log << "Writing to dec_output.sas..." << flush;
    std::ofstream output_file;
    output_file.open(file_name);
    task_dump::dump_as_SAS(*this, output_file);
    utils::g_log << "done!" << endl;
    utils::g_log << "Time for writing sas file: " << write_sas_file_timer << endl;
}

bool DecoupledRootTask::are_initial_states_consistent() const {
    for (const auto & [leaf, inner_map] : leaf_lstate_to_pvar) {
        for (const auto & [lstate, pvar] : inner_map) {
            int svar = leaf_lstate_to_svar.at(leaf).at(lstate);

            if (initial_state_values.at(svar) < initial_state_values.at(pvar))
                return false;
        }
    }
    return true;
}

vector<string> DecoupledRootTask::get_fact_names(const string &var_name) const {
    return vector<string>{
        "NegatedAtom " + var_name,
        "Atom " + var_name
    };
}

void DecoupledRootTask::create_center_variables() {
    for (int var : factoring->get_center()) {
        variables.push_back(original_root_task->variables.at(var));
        center_var_to_pvar[var] = variables.size() - 1;
    }
}

void DecoupledRootTask::create_leaf_state_variables() {
    // primary variable for leaf states
    for (int leaf = 0; leaf < factoring->get_num_leaves(); ++leaf) {
        for (int lstate = 0; lstate < factoring->get_num_leaf_states(leaf); ++lstate) {
            string name = "v(" + factoring->get_leaf_state_name(leaf, lstate) + ")";
            variables.emplace_back(name, get_fact_names(name), -1);
            leaf_lstate_to_pvar[leaf][lstate] = variables.size() - 1;
        }
    }

    // secondary variable for leaf states
    for (int leaf = 0; leaf < factoring->get_num_leaves(); ++leaf) {
        for (int lstate = 0; lstate < factoring->get_num_leaf_states(leaf); ++lstate) {
            string name = "s(" + factoring->get_leaf_state_name(leaf, lstate) + ")";
            variables.emplace_back(name, get_fact_names(name), 0);
            leaf_lstate_to_svar[leaf][lstate] = variables.size() - 1;
        }
    }
}

void DecoupledRootTask::create_goal_condition_variables() {
    // secondary variable for goal conditions
    for (const auto &g_fact : original_root_task->goals) {
        int var = g_fact.var;

        if (factoring->is_center_variable(var))
            continue;

        int leaf = factoring->get_factor(var);
        assert(leaf != -1);

        if (leaf_to_goal_svar.count(leaf) == 0) {
            string name = "g-s(" + factoring->get_leaf_name(leaf) + ")";
            variables.emplace_back(name, get_fact_names(name), 0);
            leaf_to_goal_svar[leaf] = variables.size() - 1;
        }
    }
}

// secondary variable for operator preconditions
void DecoupledRootTask::create_precondition_variables() {
    for (size_t op_id = 0; op_id < original_root_task->operators.size(); ++op_id) {
        if (!factoring->is_global_operator(op_id))
            continue;

        const auto &op = original_root_task->operators[op_id];

        // If we want to use the op name in the variable name, we need to have no white spaces!
        string no_space_op_name = op.name;
        replace(no_space_op_name.begin(), no_space_op_name.end(), ' ', '-');

        vector<vector<FactPair>> leaf_preconditions(factoring->get_num_leaves());
        for (const auto &pre : op.preconditions) {
            int var = pre.var;

            if (factoring->is_center_variable(var))
                continue;

            int leaf = factoring->get_factor(var);
            leaf_preconditions[leaf].push_back(pre);
        }

        for (size_t leaf = 0; leaf < leaf_preconditions.size(); ++leaf) {
            const vector<FactPair> &leaf_pre = leaf_preconditions[leaf];

            if (leaf_pre.empty())
                continue;

            // We have not seen this precondition and create a new secondary variable for it
            if (!same_leaf_preconditons_single_variable || precondition_to_svar.count(leaf_pre) == 0) {
                string name = "op-s(" + factoring->get_leaf_name(leaf) + "-" + no_space_op_name + ")";
                variables.emplace_back(name, get_fact_names(name), 0);
                precondition_to_svar[leaf_pre] = variables.size() - 1;
            }

            leaf_op_to_svar[leaf][op_id] = precondition_to_svar[leaf_pre];
        }
    }
}

void DecoupledRootTask::create_variables() {
    int cur_num_variables = 0;

    create_center_variables();
    utils::g_log << "\tNumber of primary center variables: " << variables.size() - cur_num_variables << endl;
    cur_num_variables = variables.size();

    create_leaf_state_variables();
    utils::g_log << "\tNumber of primary leaf state variables: " << (variables.size() - cur_num_variables) / 2 << endl;
    utils::g_log << "\tNumber of secondary leaf state variables: " << (variables.size() - cur_num_variables) / 2 << endl;
    cur_num_variables = variables.size();

    create_goal_condition_variables();
    utils::g_log << "\tNumber of secondary goal condition variables: " << variables.size() - cur_num_variables << endl;
    cur_num_variables = variables.size();

    create_precondition_variables();
    utils::g_log << "\tNumber of secondary precondition variables: " << variables.size() - cur_num_variables << endl;
    cur_num_variables = variables.size();
}

// We only keep center variables as mutexes
void DecoupledRootTask::create_mutexes() {
    mutexes.resize(variables.size());
    for (size_t var = 0; var < variables.size(); ++var)
        mutexes[var].resize(variables[var].domain_size);

    for (size_t var = 0; var < original_root_task->mutexes.size(); ++var) {
        if (!factoring->is_center_variable(var))
            continue;
        for (size_t val = 0; val < original_root_task->mutexes.at(var).size(); ++val) {
            for (const FactPair &fact : original_root_task->mutexes.at(var).at(val)) {
                if (factoring->is_center_variable(fact.var)) {
                    mutexes[center_var_to_pvar[var]][val].insert(
                        FactPair(center_var_to_pvar[fact.var], fact.value));
                }
            }
        }
    }
}

void DecoupledRootTask::create_initial_state() {
    // Fill initial state with zeros
    initial_state_values.resize(variables.size());
    fill(initial_state_values.begin(), initial_state_values.end(), 0);

    // Initial values of center variables
    for (const auto & [cvar, pvar] : center_var_to_pvar) {
        initial_state_values[pvar] = original_root_task->initial_state_values.at(cvar);
    }

    // Set the leaf state primary variables true which correspond to the initial state
    for (int leaf = 0; leaf < factoring->get_num_leaves(); ++leaf) {
        int lstate = factoring->get_initial_leaf_state(leaf);
        int pvar = leaf_lstate_to_pvar[leaf][lstate];
        initial_state_values[pvar] = 1;
    }
}

void DecoupledRootTask::create_goal() {
    vector<int> center = factoring->get_center();

    // Center variables goals
    for (const auto &g_fact : original_root_task->goals) {
        int var = g_fact.var;
        int val = g_fact.value;

        if (factoring->is_center_variable(var)) {
            int pvar = center_var_to_pvar[var];
            goals.emplace_back(pvar, val);
        }
    }

    // Leaf secondary variables
    for (const auto & [leaf, svar] : leaf_to_goal_svar) {
        goals.emplace_back(svar, 1);
    }

    // We sort the vector of goals in increasing variable order
    sort(goals.begin(), goals.end());

    assert(adjacent_find(goals.begin(), goals.end(),
                         [](const auto &a, const auto &b) {return a.var == b.var;}) == goals.end()
           && "Multiple goals for the same variable!");
}

void DecoupledRootTask::set_preconditions_of_operator(int op_id, ExplicitOperator &new_op) {
    const auto &op = original_root_task->operators[op_id];

    // Create center preconditions (map center variable ids)
    for (const auto &pre : op.preconditions) {
        int pre_var = pre.var;
        int pre_val = pre.value;

        if (factoring->is_center_variable(pre.var)) {
            int pvar = center_var_to_pvar[pre_var];
            new_op.preconditions.emplace_back(pvar, pre_val);
        }
    }

    // Add leaf preconditions (secondary variables)
    for (const auto &entry : leaf_op_to_svar) {
        int leaf = entry.first;

        // Check if we have a leaf precondition for this op
        if (leaf_op_to_svar[leaf].count(op_id) > 0) {
            int svar = leaf_op_to_svar[leaf][op_id];
            new_op.preconditions.emplace_back(svar, 1);
        }
    }

    // We sort the vector of preconditions in increasing variable order
    // sort(new_op.preconditions.begin(), new_op.preconditions.end());

    assert(op.preconditions.size() >= new_op.preconditions.size());
    assert(adjacent_find(new_op.preconditions.begin(), new_op.preconditions.end(),
                         [](const auto &a, const auto &b) {return a.var == b.var;}) == new_op.preconditions.end()
           && "Multiple preconditions for the same variable!");
}

void DecoupledRootTask::set_center_effects_of_operator(int op_id, ExplicitOperator &new_op) {
    const auto &op = original_root_task->operators[op_id];

    // Create center effects (map center variable ids)
    for (const auto &eff : op.effects) {
        assert(eff.conditions.empty());

        int eff_var = eff.fact.var;
        int eff_val = eff.fact.value;

        if (factoring->is_center_variable(eff_var)) {
            int pvar = center_var_to_pvar[eff_var];
            ExplicitEffect new_eff(pvar, eff_val, vector<FactPair>());
            new_op.effects.push_back(new_eff);
        }
    }
}

void DecoupledRootTask::set_leaf_effects_of_operator(int op_id, ExplicitOperator &new_op) {
    for (int l = 0; l < factoring->get_num_leaves(); ++l) {
        for (int ls = 0; ls < factoring->get_num_leaf_states(l); ++ls) {
            int pvar = leaf_lstate_to_pvar[l][ls];
            // cout << endl;
            // cout << factoring->get_leaf_state_name(l, ls) << " predecssors for " << op.name << ": " << endl;
            auto predecessor_ls = factoring->get_predecessors(l, ls, op_id);

            // Positive conditional effect
            for (int pred : predecessor_ls) {
                // cout << "\t" << factoring->get_leaf_state_name(l, pred) << endl;
                int svar_pred = leaf_lstate_to_svar[l][pred];
                ExplicitEffect eff(pvar, 1, vector<FactPair> {FactPair(svar_pred, 1)});
                new_op.effects.push_back(eff);
            }

            // Negative conditional effect
            ExplicitEffect eff(pvar, 0, vector<FactPair>());
            for (int pred : predecessor_ls) {
                int svar_pred = leaf_lstate_to_svar[l][pred];
                eff.conditions.emplace_back(svar_pred, 0);
            }
            // sort(eff.conditions.begin(), eff.conditions.end());

            assert(adjacent_find(eff.conditions.begin(), eff.conditions.end(),
                                 [](const auto &a, const auto &b) {return a.var == b.var;}) == eff.conditions.end()
                   && "Multiple effect conditions for the same variable!");

            new_op.effects.push_back(eff);
        }
    }
}

void DecoupledRootTask::create_operator(int op_id) {
    const auto &op = original_root_task->operators[op_id];
    assert(!op.is_an_axiom);

    ExplicitOperator new_op(op.cost, op.name, op.is_an_axiom);
    set_preconditions_of_operator(op_id, new_op);
    set_center_effects_of_operator(op_id, new_op);
    if (!implicit_effects) {
        set_leaf_effects_of_operator(op_id, new_op);
        assert(!new_op.effects.empty());
    }

    operators.push_back(new_op);
}

void DecoupledRootTask::create_operators() {
    for (size_t op_id = 0; op_id < original_root_task->operators.size(); ++op_id) {
        if (factoring->is_global_operator(op_id)) {
            create_operator(op_id);
            global_op_id_to_original_op_id[operators.size() - 1] = op_id;
            original_op_id_to_global_op_id[op_id] = operators.size() - 1;
        }
    }
    assert((int)operators.size() == factoring->get_num_global_operators());
}

void DecoupledRootTask::create_frame_axioms() {
    for (const auto & [leaf, inner_map] : leaf_lstate_to_pvar) {
        for (const auto & [lstate, pvar] : inner_map) {
            assert(leaf_lstate_to_svar.count(leaf));
            assert(leaf_lstate_to_svar[leaf].count(lstate));

            int svar = leaf_lstate_to_svar[leaf][lstate];
            string name = "ax-frame-" + factoring->get_leaf_state_name(leaf, lstate);

            ExplicitOperator new_op(0, name, true);
            new_op.preconditions.emplace_back(svar, 0);
            FactPair cond(pvar, 1);
            new_op.effects.emplace_back(ExplicitEffect(svar, 1, std::vector<FactPair> {cond}));

            axioms.push_back(new_op);
        }
    }
}

void DecoupledRootTask::create_goal_axioms() {
    for (const auto & [leaf, goal_svar] : leaf_to_goal_svar) {
        for (int goal_leaf_state : factoring->get_goal_leaf_states(leaf)) {
            string name = "ax-goal-" + factoring->get_leaf_state_name(leaf, goal_leaf_state);
            int state_svar = leaf_lstate_to_svar[leaf][goal_leaf_state];

            ExplicitOperator new_op(0, name, true);
            new_op.preconditions.emplace_back(goal_svar, 0);
            FactPair cond(state_svar, 1);
            new_op.effects.emplace_back(ExplicitEffect(goal_svar, 1, std::vector<FactPair> {cond}));

            axioms.push_back(new_op);
        }
    }
}

void DecoupledRootTask::create_precondition_axioms() {
    for (const auto & [leaf, inner_map] : leaf_op_to_svar) {
        for (const auto & [op_id, pre_svar] : inner_map) {
            assert(leaf_op_to_svar.count(leaf));
            assert(leaf_op_to_svar[leaf].count(op_id));
            assert(op_id < (int)original_root_task->operators.size());

            for (int pre_leaf_state : factoring->get_valid_precondition_leaf_states(leaf, op_id)) {
                string name = "ax-prec-" + original_root_task->operators[op_id].name + "-" +
                    factoring->get_leaf_state_name(leaf, pre_leaf_state);
                int state_svar = leaf_lstate_to_svar[leaf][pre_leaf_state];

                ExplicitOperator new_op(0, name, true);
                new_op.preconditions.emplace_back(pre_svar, 0);
                FactPair cond(state_svar, 1);
                new_op.effects.emplace_back(ExplicitEffect(pre_svar, 1, std::vector<FactPair> {cond}));

                axioms.push_back(new_op);
            }
        }
    }
}

void DecoupledRootTask::create_leaf_only_operator_axioms() {
    for (size_t op_id = 0; op_id < original_root_task->operators.size(); ++op_id) {
        if (!factoring->is_leaf_only_operator(op_id))
            continue;

        const auto &op = original_root_task->operators[op_id];
        assert(!op.effects.empty());
        int leaf = factoring->get_factor(op.effects[0].fact.var);
        assert(leaf_lstate_to_svar.count(leaf) != 0);

        // Creating center precondition
        vector<FactPair> center_conditions;
        for (const auto &pre : op.preconditions) {
            if (factoring->is_center_variable(pre.var)) {
                int center_pvar = center_var_to_pvar[pre.var];
                center_conditions.emplace_back(center_pvar, pre.value);
            }
        }

        for (const auto & [lstate, svar] : leaf_lstate_to_svar[leaf]) {
            auto predecessor_ls = factoring->get_predecessors(leaf, lstate, op_id);
            for (int pred : predecessor_ls) {
                // Trivial axiom which we can skip
                if (pred == lstate)
                    continue;

                assert(leaf_lstate_to_svar[leaf].count(pred) != 0);

                int svar_pred = leaf_lstate_to_svar[leaf][pred];
                string name = "ax-lop-(" + op.name + ")-" + factoring->get_leaf_state_name(leaf, lstate)
                    + "-" + factoring->get_leaf_state_name(leaf, pred);

                vector<FactPair> eff_conditions = center_conditions;
                eff_conditions.emplace_back(svar_pred, 1);

                ExplicitOperator new_op(0, name, true);
                new_op.preconditions.emplace_back(svar, 0);
                new_op.effects.emplace_back(svar, 1, move(eff_conditions));

                axioms.push_back(new_op);
            }
        }
    }
}

void DecoupledRootTask::create_axioms() {
    create_frame_axioms();
    assert((int)axioms.size() == factoring->get_num_all_leaf_states());

    create_goal_axioms();
    assert((int)axioms.size() == factoring->get_num_all_leaf_states() +
           factoring->get_num_all_goal_leaf_states());

    create_precondition_axioms();

    create_leaf_only_operator_axioms();

    assert(ranges::all_of(axioms, [](const auto &axiom)
                          {return axiom.effects.size() == 1;}));
    assert(ranges::all_of(axioms, [](const auto &axiom)
                          {return axiom.preconditions.size() == 1;}));
    assert(ranges::all_of(axioms, [](const auto &axiom)
                          {return axiom.preconditions.at(0).var == axiom.effects.at(0).fact.var;}));
    assert(ranges::all_of(axioms, [](const auto &axiom)
                          {return axiom.preconditions.at(0).value != axiom.effects.at(0).fact.value;}));
}

// TODO: release more memory
void DecoupledRootTask::release_memory() {
    precondition_to_svar.clear();
}

void DecoupledRootTask::dump() const {
    task_properties::dump_task(TaskProxy(*this), true, true);
}


shared_ptr<AbstractTask> DecoupledRootTask::get_original_root_task() const {
    return original_root_task;
}

void DecoupledRootTask::set_center_values(const State &dec_state, vector<int> &state) const {
    for (auto const & [original_var, decoupled_pvar] : center_var_to_pvar) {
        state[original_var] = dec_state[decoupled_pvar].get_value();
    }
}

void DecoupledRootTask::set_random_leave_values(const State &dec_state, vector<int> &state) const {
    for (int l = 0; l < factoring->get_num_leaves(); ++l) {
        vector<int> reached_leaf_states;
        for (const auto & [lstate, svar] : leaf_lstate_to_svar.at(l)) {
            if (dec_state[svar].get_value() == 1) {
                reached_leaf_states.push_back(lstate);
            }
        }
        assert(!reached_leaf_states.empty());
        int selected_leaf_state = *rng->choose(reached_leaf_states);
        factoring->add_leaf_facts_to_state(state, l, selected_leaf_state);
    }
}

bool DecoupledRootTask::is_valid_decoupled_state(const State &dec_state) const {
    for (int l = 0; l < factoring->get_num_leaves(); ++l) {
        vector<int> reached_leaf_states;
        for (const auto & [lstate, svar] : leaf_lstate_to_svar.at(l)) {
            if (dec_state[svar].get_value() == 1) {
                reached_leaf_states.push_back(lstate);
            }
        }

        if (reached_leaf_states.empty()) {
            return false;
        }
    }
    return true;
}

const ExplicitEffect &DecoupledRootTask::get_effect(int op_id, int effect_id, bool is_axiom) const {
    if (!implicit_effects || is_axiom)
        return RootTask::get_effect(op_id, effect_id, is_axiom);

    if (effect_id < RootTask::get_num_operator_effects(op_id, is_axiom))
        return RootTask::get_effect(op_id, effect_id, is_axiom);
}

int DecoupledRootTask::get_num_operator_effects(int op_index, bool is_axiom) const {
    if (!implicit_effects || is_axiom)
        return RootTask::get_num_operator_effects(op_index, is_axiom);
}

int DecoupledRootTask::get_num_operator_effect_conditions(int op_index, int eff_index, bool is_axiom) const {
    if (!implicit_effects || is_axiom)
        return RootTask::get_num_operator_effect_conditions(op_index, eff_index, is_axiom);

    if (eff_index < RootTask::get_num_operator_effects(op_index, is_axiom))
        return RootTask::get_num_operator_effect_conditions(op_index, eff_index, is_axiom);
}

FactPair DecoupledRootTask::get_operator_effect_condition(int op_index, int eff_index, int cond_index, bool is_axiom) const {
    if (!implicit_effects || is_axiom)
        return RootTask::get_operator_effect_condition(op_index, eff_index, cond_index, is_axiom);

    if (eff_index < RootTask::get_num_operator_effects(op_index, is_axiom))
        return RootTask::get_operator_effect_condition(op_index, eff_index, cond_index, is_axiom);
}

FactPair DecoupledRootTask::get_operator_effect(int op_index, int eff_index, bool is_axiom) const {
    if (!implicit_effects || is_axiom)
        return RootTask::get_operator_effect(op_index, eff_index, is_axiom);

    if (eff_index < RootTask::get_num_operator_effects(op_index, is_axiom))
        return RootTask::get_operator_effect(op_index, eff_index, is_axiom);
}

class DecoupledRootTaskFeature : public plugins::TypedFeature<AbstractTask, DecoupledRootTask> {
public:
    DecoupledRootTaskFeature() : TypedFeature("decoupled") {
        document_title("Decoupled task");
        document_synopsis(
            "A decoupled transformation of the root task.");
        utils::add_rng_options(*this);

        add_option<shared_ptr<decoupling::Factoring>>("factoring",
                                                      "method that computes the factoring");
        add_option<bool>("same_leaf_preconditons_single_variable", "The same preconditions of leaf have a single secondary variables.", "true");
        add_option<bool>("implicit_effects", "Represent effects implicitly and create them on demand", "false");
        add_option<bool>("write_sas_file", "Writes the decoupled task to dec_output.sas and terminates.", "false");
        add_option<bool>("dump_task", "Dumps the task to the console", "false");
    }

    virtual shared_ptr<DecoupledRootTask> create_component(const plugins::Options &options, const utils::Context &) const override {
        return make_shared<DecoupledRootTask>(options);
    }
};

static plugins::FeaturePlugin<DecoupledRootTaskFeature> _plugin;
}

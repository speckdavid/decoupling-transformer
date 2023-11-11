#include "decoupled_root_task.h"

#include "../axioms.h"
#include "../decoupling/factoring.h"
#include "../plugins/plugin.h"
#include "../task_utils/task_dump.h"
#include "../task_utils/task_properties.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>

using namespace std;

namespace tasks {
DecoupledRootTask::DecoupledRootTask(const plugins::Options &options)
    : RootTask(),
      original_root_task(dynamic_pointer_cast<RootTask>(tasks::g_root_task)),
      factoring(options.get<shared_ptr<decoupling::Factoring>>("factoring")) {
    factoring->compute_factoring();

    utils::Timer transformation_timer;

    // TODO: statistics: time, number axioms, number variables, etc.
    // TODO: options to only dump to sas

    TaskProxy original_task_proxy(*original_root_task);
    task_properties::verify_no_axioms(original_task_proxy);
    task_properties::verify_no_conditional_effects(original_task_proxy);

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

    // task_properties::dump_task(task_proxy);

    utils::g_log << "Time for decoupled transformation: " << transformation_timer << endl;
    print_statistics();

    if (options.get<bool>("write_sas_file")) {
        write_sas_file("dec_output.sas");
        utils::exit_with(utils::ExitCode::SEARCH_UNSOLVED_INCOMPLETE);
    }
}

void DecoupledRootTask::print_statistics() const {
    int num_primary_vars = count_if(variables.begin(), variables.end(), [](const auto &var)
                                    {return var.axiom_layer == -1;});
    int num_secondary_vars = variables.size() - num_primary_vars;

    utils::g_log << "Original task size: " << original_root_task->get_encoding_size(false) << endl;
    utils::g_log << "Original number of variables: " << original_root_task->variables.size() << endl;
    utils::g_log << "Original number of primary variables: " << original_root_task->variables.size() << endl;
    utils::g_log << "Original number of secondary variables: " << 0 << endl;
    utils::g_log << "Original number of operators: " << original_root_task->operators.size() << endl;
    utils::g_log << "Original number of axioms: " << original_root_task->axioms.size() << endl;

    utils::g_log << "Task size: " << get_encoding_size(false) << endl;
    utils::g_log << "Number of variables: " << variables.size() << endl;
    utils::g_log << "Number of primary variables: " << num_primary_vars << endl;
    utils::g_log << "Number of secondary variables: " << num_secondary_vars << endl;
    utils::g_log << "Number of operators: " << operators.size() << endl;
    utils::g_log << "Number of axioms: " << axioms.size() << endl;
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

void DecoupledRootTask::create_variables() {
    vector<int> center = factoring->get_center();
    for (int var : center) {
        variables.push_back(original_root_task->variables.at(var));
        center_var_to_pvar[var] = variables.size() - 1;
    }

    // primary variable for leaf states
    for (int leaf = 0; leaf < factoring->get_num_leaves(); ++leaf) {
        for (int lstate = 0; lstate < factoring->get_num_leaf_states(leaf); ++lstate) {
            string name = "v(" + factoring->get_leaf_state_name(leaf, lstate) + ")";
            variables.emplace_back(name, vector<string>{"False", "True"}, -1);
            leaf_lstate_to_pvar[leaf][lstate] = variables.size() - 1;
        }
    }

    // secondary variable for leaf states
    for (int leaf = 0; leaf < factoring->get_num_leaves(); ++leaf) {
        for (int lstate = 0; lstate < factoring->get_num_leaf_states(leaf); ++lstate) {
            string name = "s(" + factoring->get_leaf_state_name(leaf, lstate) + ")";
            variables.emplace_back(name, vector<string>{"False", "True"}, 0);
            leaf_lstate_to_svar[leaf][lstate] = variables.size() - 1;
        }
    }

    // secondary variable for goal conditions
    for (const auto &g_fact : original_root_task->goals) {
        int var = g_fact.var;

        if (factoring->is_center_variable(var))
            continue;

        int leaf = factoring->get_factor(var);
        assert(leaf != -1);

        if (leaf_to_goal_svar.count(leaf) == 0) {
            string name = "g-s(" + factoring->get_leaf_name(leaf) + ")";
            variables.emplace_back(name, vector<string>{"False", "True"}, 0);
            leaf_to_goal_svar[leaf] = variables.size() - 1;
        }
    }

    // secondary variable for operator preconditions
    int op_id = 0;
    for (const auto &op : original_root_task->operators) {
        if (!factoring->is_global_operator(op_id))
           continue;

        // If we want to use the op name in the variable name, we need to have no white spaces!
        string no_space_op_name = op.name;
        replace(no_space_op_name.begin(), no_space_op_name.end(), ' ', '-');

        for (const auto &pre : op.preconditions) {
            int var = pre.var;


            if (factoring->is_center_variable(var))
                continue;

            int leaf = factoring->get_factor(var);
            assert(leaf != -1);

            if (leaf_op_to_svar[leaf].count(op_id) == 0) {
                string name = "op-s(" + factoring->get_leaf_name(leaf) + "-" + no_space_op_name + ")";
                variables.emplace_back(name, vector<string>{"False", "True"}, 0);
                leaf_op_to_svar[leaf][op_id] = variables.size() - 1;
            }
        }
        ++op_id;
    }
}

// TODO: At the moment we simply ignore mutexes by leaving them empty.
// At least we could keep the mutexes over the center variables
void DecoupledRootTask::create_mutexes() {
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

void DecoupledRootTask::set_precondition_of_operator(int op_id, ExplicitOperator &new_op) {
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
    sort(new_op.preconditions.begin(), new_op.preconditions.end());

    assert(op.preconditions.size() >= new_op.preconditions.size());
    assert(adjacent_find(new_op.preconditions.begin(), new_op.preconditions.end(),
                         [](const auto &a, const auto &b) {return a.var == b.var;}) == new_op.preconditions.end()
           && "Multiple preconditions for the same variable!");
}

void DecoupledRootTask::set_effect_of_operator(int op_id, ExplicitOperator &new_op) {
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

    // Effects on leaf states
    for (int l = 0; l < factoring->get_num_leaves(); ++l) {
        for (int ls = 0; ls < factoring->get_num_leaf_states(l); ++ls) {
            int pvar = leaf_lstate_to_pvar[l][ls];
            auto predecessor_ls = factoring->get_predecessors(l, ls, op_id);

            // Positive conditional effect
            for (int pred : predecessor_ls) {
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
            sort(eff.conditions.begin(), eff.conditions.end());

            assert(adjacent_find(eff.conditions.begin(), eff.conditions.end(),
                                 [](const auto &a, const auto &b) {return a.var == b.var;}) == eff.conditions.end()
                   && "Multiple effect conditions for the same variable!");

            new_op.effects.push_back(eff);
        }

        // Sort the vector based on the 'eff'. Is this necessary?
        sort(new_op.effects.begin(), new_op.effects.end(), [](const ExplicitEffect &a, const ExplicitEffect &b) {
                 return a.fact < b.fact;
             });
    }

    assert(!new_op.effects.empty());
}

void DecoupledRootTask::create_operator(int op_id) {
    const auto &op = original_root_task->operators[op_id];
    assert(!op.is_an_axiom);

    ExplicitOperator new_op(op.cost, op.name, op.is_an_axiom);
    set_precondition_of_operator(op_id, new_op);
    set_effect_of_operator(op_id, new_op);

    operators.push_back(new_op);
}

void DecoupledRootTask::create_operators() {
    for (size_t op_id = 0; op_id < original_root_task->operators.size(); ++op_id) {
        if (factoring->is_global_operator(op_id)) {
            create_operator(op_id);
        }
    }

    // for (size_t i = 0; i < operators.size(); ++i) {
    //     utils::g_log << operators[i].name << ":" << endl;
    //     utils::g_log << "\tpre: " << operators[i].preconditions << endl;
    //     utils::g_log << "\teff: " << endl;
    //     for (auto const &eff : operators[i].effects) {
    //         utils::g_log << "\t   " << eff.fact << " if " << eff.conditions << endl;
    //     }
    //     cout << endl;
    // }

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

    // for (const ExplicitOperator &axiom : axioms) {
    //     vector<FactPair> conds = axiom.effects.at(0).conditions;
    //     conds.insert(conds.begin(), axiom.preconditions.at(0));
    //     utils::g_log << axiom.name << ": ("
    //                  << axiom.effects.at(0).fact << ") <= "
    //                  << conds << endl;
    // }
}

class DecoupledRootTaskFeature : public plugins::TypedFeature<AbstractTask, DecoupledRootTask> {
public:
    DecoupledRootTaskFeature() : TypedFeature("decoupled") {
        document_title("Decoupled task");
        document_synopsis(
            "A decoupled transformation of the root task.");

        add_option<shared_ptr<decoupling::Factoring>>("factoring",
                                                      "method that computes the factoring");
        add_option<bool>("write_sas_file", "Writes the decoupled task to dec_output.sas and terminates.", "false");
    }

    virtual shared_ptr<DecoupledRootTask> create_component(const plugins::Options &options, const utils::Context &) const override {
        return make_shared<DecoupledRootTask>(options);
    }
};

static plugins::FeaturePlugin<DecoupledRootTaskFeature> _plugin;
}

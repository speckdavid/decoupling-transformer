#include "decoupled_root_task.h"

#include "../axioms.h"
#include "../decoupling/factoring.h"
#include "../plugins/plugin.h"
#include "../task_utils/task_dump.h"
#include "../task_utils/task_properties.h"

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

    TaskProxy original_task_proxy(*original_root_task);
    task_properties::verify_no_axioms(original_task_proxy);
    task_properties::verify_no_conditional_effects(original_task_proxy);

    create_variables();
    create_mutexes();
    create_operators();
    create_initial_state();
    create_goal();
    create_axioms();

    TaskProxy task_proxy(*this);

    // This is also done in the root task which is honestly quite hacky!
    AxiomEvaluator &axiom_evaluator = g_axiom_evaluators[task_proxy];
    axiom_evaluator.evaluate(initial_state_values);

    // task_properties::dump_task(task_proxy);
    exit(0);
}

bool DecoupledRootTask::are_initial_states_consistent(bool exact_match) const {
    vector<int> initial_state_values_copy = initial_state_values;
    AxiomEvaluator &axiom_evaluator = g_axiom_evaluators[TaskProxy(*this)];
    axiom_evaluator.evaluate(initial_state_values_copy);

    for (const auto & [leaf, inner_map] : leaf_lstate_to_pvar) {
        for (const auto & [lstate, pvar] : inner_map) {
            int svar = leaf_lstate_to_svar.at(leaf).at(lstate);

            if (exact_match) {
                if (initial_state_values_copy.at(pvar) != initial_state_values_copy.at(svar))
                    return false;
            } else {
                if (initial_state_values_copy.at(pvar) < initial_state_values_copy.at(svar))
                    return false;
            }
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
            string name = "v(" + to_string(leaf) + "," + to_string(lstate) + ")";
            variables.emplace_back(name, vector<string>{"False", "True"}, -1);
            leaf_lstate_to_pvar[leaf][lstate] = variables.size() - 1;
        }
    }

    // secondary variable for leaf states
    for (int leaf = 0; leaf < factoring->get_num_leaves(); ++leaf) {
        for (int lstate = 0; lstate < factoring->get_num_leaf_states(leaf); ++lstate) {
            string name = "s(" + to_string(leaf) + "," + to_string(lstate) + ")";
            variables.emplace_back(name, vector<string>{"False", "True"}, 0);
            leaf_lstate_to_svar[leaf][lstate] = variables.size() - 1;
        }
    }

    // secondary variable for goal conditions
    for (const auto &g_fact : original_root_task->goals) {
        int var = g_fact.var;

        if (factoring->is_center_variable(var))
            continue;

        int leaf = factoring->get_leaf_of_variable(var);
        assert(leaf != -1);

        if (leaf_to_goal_svar.count(leaf) == 0) {
            string name = "g-s(" + to_string(leaf) + ")";
            variables.emplace_back(name, vector<string>{"False", "True"}, 0);
            leaf_to_goal_svar[leaf] = variables.size() - 1;
        }
    }

    // secondary variable for operator preconditions
    int op_id = 0;
    for (const auto &op : original_root_task->operators) {
        for (const auto &pre : op.preconditions) {
            int var = pre.var;

            if (factoring->is_center_variable(var))
                continue;

            int leaf = factoring->get_leaf_of_variable(var);
            assert(leaf != -1);

            if (leaf_op_to_svar[leaf].count(op_id) == 0) {
                string name = "op-s(" + to_string(leaf) + "," + to_string(op_id) + ")";
                variables.emplace_back(name, vector<string>{"False", "True"}, 0);
                leaf_op_to_svar[leaf][op_id] = variables.size() - 1;
            }

            // Collect the preconditions for each leaf and op
            leaf_op_to_pre[leaf][op_id].push_back(pre);
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

    // Copy center preconditions
    for (const auto &pre : op.preconditions) {
        if (factoring->is_center_variable(pre.var))
            new_op.preconditions.push_back(pre);
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

    // Copy center effects
    for (const auto &eff : op.effects) {
        assert(eff.conditions.empty());

        if (factoring->is_center_variable(eff.fact.var))
            new_op.effects.push_back(eff);
    }

    // Effects on leaf states
    for (int l = 0; l < factoring->get_num_leaves(); ++l) {
        for (int ls = 0; ls < factoring->get_num_leaf_states(l); ++ls) {
            int pvar = leaf_lstate_to_pvar[l][ls];
            set<int> predecessor_ls = factoring->get_predecessors(l, ls, op_id);

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
        create_operator(op_id);
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

    assert(operators.size() == original_root_task->operators.size());
}

void DecoupledRootTask::create_frame_axioms() {
    for (const auto & [leaf, inner_map] : leaf_lstate_to_pvar) {
        for (const auto & [lstate, pvar] : inner_map) {
            assert(leaf_lstate_to_svar.count(leaf));
            assert(leaf_lstate_to_svar[leaf].count(lstate));

            int svar = leaf_lstate_to_svar[leaf][lstate];
            string name = "frame-" + to_string(leaf) + "-" + to_string(lstate);

            ExplicitOperator new_op(0, name, true);
            FactPair cond(pvar, 1);
            new_op.effects.emplace_back(ExplicitEffect(svar, 1, std::vector<FactPair> {cond}));

            axioms.push_back(new_op);
        }
    }
}

void DecoupledRootTask::create_precondition_axioms() {
    // TODO: implement
}


void DecoupledRootTask::create_goal_axioms() {
    // TODO: implement
}

void DecoupledRootTask::create_leaf_only_operator_axioms() {
    for (size_t op_id = 0; op_id < original_root_task->operators.size(); ++op_id) {
        if (!factoring->is_leaf_only_operator(op_id))
            continue;

        const auto &op = original_root_task->operators[op_id];
        assert(!op.effects.empty());
        int leaf = factoring->get_leaf_of_variable(op.effects[0].fact.var);

        // Creating center precondition
        vector<FactPair> center_conditions;
        for (const auto &pre : op.preconditions) {
            if (factoring->is_center_variable(pre.var))
                center_conditions.push_back(pre);
        }

        for (const auto & [lstate, svar] : leaf_lstate_to_svar[leaf]) {
            set<int> predecessor_ls = factoring->get_predecessors(leaf, lstate, op_id);
            for (int pred : predecessor_ls) {
                int svar_pred = leaf_lstate_to_svar[lstate][pred];
                string name = "lop-" + op.name + "-" + to_string(leaf) + "-" + to_string(lstate) + "-" + to_string(pred);

                vector<FactPair> eff_conditions = center_conditions;
                eff_conditions.emplace_back(svar_pred, 1);
                assert(eff_conditions.size() == center_conditions.size() + 1);

                ExplicitOperator new_op(0, name, true);
                new_op.effects.emplace_back(svar, 1, move(eff_conditions));

                axioms.push_back(new_op);
            }
        }
    }
}

void DecoupledRootTask::create_axioms() {
    create_frame_axioms();

    // The exact match is important to be done before the leaf only operator axioms are created
    assert(are_initial_states_consistent(true));

    create_leaf_only_operator_axioms();

    assert(are_initial_states_consistent(false));
}

class DecoupledRootTaskFeature : public plugins::TypedFeature<AbstractTask, DecoupledRootTask> {
public:
    DecoupledRootTaskFeature() : TypedFeature("decoupled") {
        document_title("Decoupled task");
        document_synopsis(
            "A decoupled transformation of the root task.");

        add_option<shared_ptr<decoupling::Factoring>>("factoring",
                                                      "method that computes the factoring");
    }

    virtual shared_ptr<DecoupledRootTask> create_component(const plugins::Options &options, const utils::Context &) const override {
        return make_shared<DecoupledRootTask>(options);
    }
};

static plugins::FeaturePlugin<DecoupledRootTaskFeature> _plugin;
}

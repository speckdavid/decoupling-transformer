#include "decoupled_root_task.h"

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

    TaskProxy task_proxy(*original_root_task);
    task_properties::verify_no_axioms(task_proxy);
    task_properties::verify_no_conditional_effects(task_proxy);

    create_variables();
    create_mutexes();
    create_operators();
    create_initial_state();
    create_goal();
    create_operators();
    create_axioms();

    // task_properties::dump_task(TaskProxy(*this));
    exit(0);
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
            string name = "d(" + to_string(leaf) + "," + to_string(lstate) + ")";
            variables.emplace_back(name, vector<string>{"False", "True"}, 0);
            leaf_lstate_to_dvar[leaf][lstate] = variables.size() - 1;
        }
    }

    // secondary variable for goal conditions
    for (const auto &g_fact : original_root_task->goals) {
        int var = g_fact.var;

        if (factoring->is_center_variable(var))
            continue;

        int leaf = factoring->get_leaf_of_variables(var);
        assert(leaf != -1);

        if (leaf_to_goal_pvar.count(leaf) == 0) {
            string name = "g-d(" + to_string(leaf) + ")";
            variables.emplace_back(name, vector<string>{"False", "True"}, 0);
            leaf_to_goal_pvar[leaf] = variables.size() - 1;
        }
    }

    // TODO: create derived variables for operator preconditions!
}

// TODO: At the moment we simply ignore mutexes by leaving them empty.
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
    for (const auto & [leaf, pvar] : leaf_to_goal_pvar) {
        goals.emplace_back(pvar, 1);
    }

    // We sort the vector of goals in increasing variable order
    sort(goals.begin(), goals.end());

    assert(adjacent_find(goals.begin(), goals.end(),
                         [](const auto &a, const auto &b) {return a.var == b.var;}) == goals.end()
           && "Multiple goals for the same variable!");
}

void DecoupledRootTask::create_operators() {
}

void DecoupledRootTask::create_axioms() {
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

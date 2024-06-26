#include "permutation_group_task.h"

#include "../structural_symmetries/group.h"
#include "../structural_symmetries/permutation.h"
#include "../utils/logging.h"
#include "../utils/system.h"
#include "../task_utils/causal_graph.h"
#include "../task_utils/task_properties.h"

using namespace std;
using namespace structural_symmetries;

namespace extra_tasks {
PermutationGroupTask::PermutationGroupTask(
    const shared_ptr<AbstractTask> &parent,
    const shared_ptr<Group> &group)
    : DelegatingTask(parent),
      group(group) {
    TaskProxy parent_task_proxy(*parent);
    task_properties::verify_no_axioms(parent_task_proxy);
    task_properties::verify_no_conditional_effects(parent_task_proxy);

    if (!group->is_initialized()) {
        group->compute_symmetries(TaskProxy(*parent));
        if (!group->has_symmetries()) {
            utils::g_log << "No symmetries found, aborting.." << endl;
            utils::exit_with(utils::ExitCode::SEARCH_UNSOLVED_INCOMPLETE);
        }
    }


    for (int g = 0; g < group->get_num_generators(); ++g) {
        vector<tasks::ExplicitEffect> effects;
        auto const &gen = group->get_permutation(g);
        for (int var : gen.get_affected_vars()) {
            for (int val = 0; val < get_variable_domain_size(var); ++val) {
                auto [new_var, new_val] = gen.get_new_var_val_by_old_var_val(var, val);
                effects.emplace_back(new_var, new_val, vector<FactPair>{FactPair(var, val)});

                effect_variables.insert(new_var);
                condition_variables.insert(var);
            }
        }
        generator_effects.push_back(effects);
    }

    relevant_variables = effect_variables;
    relevant_variables.insert(condition_variables.begin(), condition_variables.end());

    set<int> open_vars = relevant_variables;
    while (!open_vars.empty()) {
        int var = *open_vars.begin();
        open_vars.erase(var);
        set<int> component;
        compute_causal_graph_component_of_variable(var, component);
        relevant_components.push_back(component);
        for (int var : component) {
            open_vars.erase(var);
        }
    }

    auto pg_cg = causal_graph::get_causal_graph(this);
    for (int var = 0; var < get_num_variables(); ++var) {
        set<int> neighbors(pg_cg.get_predecessors(var).begin(), pg_cg.get_predecessors(var).end());
        neighbors.insert(pg_cg.get_successors(var).begin(), pg_cg.get_successors(var).end());
        variable_degree.push_back(neighbors.size());
    }
}

void PermutationGroupTask::compute_causal_graph_component_of_variable(int start_var, set<int> &component) const {
    // Ensure the component set is empty and the start variable is relevant
    assert(component.empty());
    assert(relevant_variables.count(start_var));

    // Retrieve the causal graph associated with this task
    auto pg_cg = causal_graph::get_causal_graph(this);

    // Initialize a set for variables to explore
    set<int> open;
    open.insert(start_var);

    // Explore the causal graph using a breadth-first search
    while (!open.empty()) {
        // Get the next variable and mark it as processed
        int var = *open.begin();
        open.erase(var);
        component.insert(var);

        // Check successors and add to the open set if relevant and not already processed
        for (int succ : pg_cg.get_successors(var)) {
            if (component.count(succ) == 0 && relevant_variables.count(succ) > 0) {
                open.insert(succ);
            }
        }

        // Check predecessors and add to the open set if relevant and not already processed
        for (int pred : pg_cg.get_predecessors(var)) {
            if (component.count(pred) == 0 && relevant_variables.count(pred) > 0) {
                open.insert(pred);
            }
        }
    }
}

int PermutationGroupTask::get_operator_cost(int /*index*/, bool /*is_axiom*/) const {
    return 1;
}

string PermutationGroupTask::get_operator_name(int index, bool /*is_axiom*/) const {
    return "Generator " + to_string(index);
}

int PermutationGroupTask::get_num_operators() const {
    return group->get_num_generators();
}

int PermutationGroupTask::get_num_operator_preconditions(int /*index*/, bool /*is_axiom*/) const {
    return 0;
}

FactPair PermutationGroupTask::get_operator_precondition(int /*op_index*/, int /*fact_index*/, bool /*is_axiom*/) const {
    utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    return FactPair::no_fact;
}

int PermutationGroupTask::get_num_operator_effects(int op_index, bool /*is_axiom*/) const {
    return generator_effects.at(op_index).size();
}

int PermutationGroupTask::get_num_operator_effect_conditions(
    int op_index, int eff_index, bool /*is_axiom*/) const {
    return generator_effects.at(op_index).at(eff_index).conditions.size();
}

FactPair PermutationGroupTask::get_operator_effect_condition(
    int op_index, int eff_index, int cond_index, bool /*is_axiom*/) const {
    return generator_effects.at(op_index).at(eff_index).conditions.at(cond_index);
}

FactPair PermutationGroupTask::get_operator_effect(
    int op_index, int eff_index, bool /*is_axiom*/) const {
    return generator_effects.at(op_index).at(eff_index).fact;
}

int PermutationGroupTask::convert_operator_index_to_parent(int /*index*/) {
    return OperatorID::no_operator.get_index();
}

void PermutationGroupTask::write_causal_graph(const string &file_name) const {
    auto pg_cg = causal_graph::get_causal_graph(this);
    pg_cg.to_dot(TaskProxy(*this), file_name);
}

const set<int> &PermutationGroupTask::get_effect_variables() const {
    return effect_variables;
}
const set<int> &PermutationGroupTask::get_condition_variables() const {
    return condition_variables;
}
const set<int> &PermutationGroupTask::get_relevant_variables() const {
    return relevant_variables;
}

void PermutationGroupTask::get_relevant_components(std::vector<std::vector<int>> &components,
                                                   bool components_sorted_by_degree) const {
    assert(components.empty());

    auto pg_cg = causal_graph::get_causal_graph(this);

    for (const auto &component : relevant_components) {
        std::vector<int> cur(component.begin(), component.end());
        if (components_sorted_by_degree) {
            sort(cur.begin(), cur.end(), [&](int a, int b) {
                     int degree_a = variable_degree[a];
                     int degree_b = variable_degree[b];

                     if (degree_a != degree_b) {
                         return degree_a > degree_b; // Sort by degree in decreasing order
                     } else {
                         return a < b; // Sort by ID in increasing order for ties
                     }
                 });
        }
        components.push_back(cur);
    }
}
}

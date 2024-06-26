#ifndef TASKS_PERMUTATION_GROUP_TASK_H
#define TASKS_PERMUTATION_GROUP_TASK_H

#include "delegating_task.h"

#include "root_task.h"

#include <set>
#include <vector>

namespace structural_symmetries {
class Group;
class Permutation;
}

namespace extra_tasks {
class PermutationGroupTask : public tasks::DelegatingTask {
protected:
    const std::shared_ptr<structural_symmetries::Group> group;
    std::vector<std::vector<tasks::ExplicitEffect>> generator_effects;
    std::set<int> condition_variables;
    std::set<int> effect_variables;
    std::set<int> relevant_variables;

    std::vector<std::set<int>> relevant_components;
    std::vector<int> variable_degree;

    void compute_causal_graph_component_of_variable(int start_var, std::set<int> &component) const;


public:
    PermutationGroupTask(
        const std::shared_ptr<AbstractTask> &parent,
        const std::shared_ptr<structural_symmetries::Group> &group);
    ~PermutationGroupTask() = default;

    virtual int get_operator_cost(int index, bool is_axiom) const override;
    virtual std::string get_operator_name(int index, bool is_axiom) const override;
    virtual int get_num_operators() const override;
    virtual int get_num_operator_preconditions(int index, bool is_axiom) const override;
    virtual FactPair get_operator_precondition(
        int op_index, int fact_index, bool is_axiom) const override;
    virtual int get_num_operator_effects(int op_index, bool is_axiom) const override;
    virtual int get_num_operator_effect_conditions(
        int op_index, int eff_index, bool is_axiom) const override;
    virtual FactPair get_operator_effect_condition(
        int op_index, int eff_index, int cond_index, bool is_axiom) const override;
    virtual FactPair get_operator_effect(
        int op_index, int eff_index, bool is_axiom) const override;
    virtual int convert_operator_index_to_parent(int index);

    // Additional functionality
    void write_causal_graph(const std::string &file_name) const;

    const std::set<int> &get_effect_variables() const;
    const std::set<int> &get_condition_variables() const;
    const std::set<int> &get_relevant_variables() const;

    // Computes variable components based in the causal graph.
    // Note: 1) We treat edges as undirected
    //       2) We only considers relevant variables (= variables part of any generator)
    void get_relevant_components(std::vector<std::vector<int>> &components,
                                 bool components_sorted_by_degree) const;
};
}

#endif

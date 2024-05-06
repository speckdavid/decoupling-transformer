#ifndef TASKS_PERMUTATION_GROUP_TASK_H
#define TASKS_PERMUTATION_GROUP_TASK_H

#include "delegating_task.h"

#include "root_task.h"

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
};
}

#endif

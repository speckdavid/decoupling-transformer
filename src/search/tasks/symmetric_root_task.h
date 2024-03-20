#ifndef TASKS_SYMMETRIC_ROOT_TASK_H
#define TASKS_SYMMETRIC_ROOT_TASK_H

#include "root_task.h"

#include "../task_proxy.h"

#include <map>
#include <unordered_map>

namespace plugins {
class Options;
}

namespace structural_symmetries {
class Group;
class Permutation;
}

namespace tasks {

enum EmptyValueStrategy {NONE = 0, INIT = 1, RANDOM = 2};

/*
  Task transformation that encodes symmetry pruning into condition effects
*/
class SymmetricRootTask : public RootTask {
    std::shared_ptr<RootTask> original_root_task; // TODO check if this is needed
    std::shared_ptr<structural_symmetries::Group> group;
    EmptyValueStrategy empty_value_strategy;
    bool compute_perfect_canonical;

    std::unique_ptr<structural_symmetries::Permutation> initial_state_permutation;
    std::vector<int> base_state_for_op_permutation;

    std::vector<int> get_state_for_operator_permutation(const ExplicitOperator &op) const;

public:
    SymmetricRootTask(const plugins::Options &options);
    virtual ~SymmetricRootTask() override = default;

    void reconstruct_plan_if_necessary(std::vector<OperatorID> &path,
                                       std::vector<State> &states,
                                       StateRegistry &registry) const override;

    virtual TaskProxy get_task_proxy_for_plan_saving() const {
        // TODO check if this is indeed needed
        // If we run symmetry search, we need the original task to save the reconstructed plan.
        return TaskProxy(*original_root_task);
    }

    // TODO check if this is needed
    std::shared_ptr<AbstractTask> get_original_root_task() const;

protected:
    void print_statistics() const;
    void write_sas_file(const std::string &file_name) const;

    void create_initial_state();

    // operators
    void set_symmetry_effects_of_operator(int op_id, ExplicitOperator &op);
    void create_operator(int op_id);
    void create_operators();

    void release_memory();

    void dump() const;
};
}

#endif

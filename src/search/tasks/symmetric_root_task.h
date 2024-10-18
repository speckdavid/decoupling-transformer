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

enum EmptyValueStrategy {NONE = 0, INIT = 1, RANDOM = 2, GOAL = 3, INIT_GOAL = 4, SPLIT_CONTEXT = 5};

/*
  Task transformation that encodes symmetry pruning into condition effects
*/
class SymmetricRootTask : public RootTask {
    std::shared_ptr<RootTask> original_root_task;
    std::shared_ptr<structural_symmetries::Group> group;
    EmptyValueStrategy empty_value_strategy;
    bool skip_mutex_preconditions;
    bool skip_unaffected_variables;
    bool skip_unaffected_variables_relevant_permutations;
    bool decoupled_splitting;
    int max_number_contexts_per_operator;
    bool compute_perfect_canonical;

    std::unique_ptr<structural_symmetries::Permutation> initial_state_permutation;
    std::vector<int> base_state_for_op_permutation;
    std::vector<int> new_op_id_to_original_op_id; // only used for empty_value_strategy==SPLIT_CONTEXT and decoupled_splitting=false
    std::vector<std::vector<int>> decoupled_splitting_implied_relevant_vars;

    std::vector<int> get_operator_post_condition(const ExplicitOperator &op) const;

    void compute_decoupled_splitting_implied_relevant_vars();

    std::vector<int> get_split_variables(const ExplicitOperator &op) const;

public:
    explicit SymmetricRootTask(const plugins::Options &options);
    
    ~SymmetricRootTask() override = default;

    void reconstruct_plan_if_necessary(std::vector<OperatorID> &path,
                                       std::vector<State> &states,
                                       StateRegistry &registry) const override;

    virtual TaskProxy get_task_proxy_for_plan_saving() const override {
        // TODO check if this is indeed needed
        // If we run symmetry search, we need the original task to save the reconstructed plan.
        return TaskProxy(*original_root_task);
    }

protected:
    void print_statistics() const;
    void write_sas_file(const std::string &file_name) const;

    void create_initial_state();

    // operators
    std::unique_ptr<structural_symmetries::Permutation> get_permutation_for_operator(
            const ExplicitOperator &op) const;

    void set_symmetry_effects_of_operator(
            const ExplicitOperator &orig_op,
            ExplicitOperator &new_op,
            const std::unique_ptr<structural_symmetries::Permutation> &perm) const;

    void add_context_split_cond_effs_recursive(
            size_t var_id,
            std::vector<FactPair> &cond_eff_preconditions,
            const std::vector<int> &component_split_vars,
            const std::vector<ExplicitEffect> &component_effects,
            ExplicitOperator &new_op);

    void create_operators_context_split_recursive(
            size_t var_id,
            std::vector<FactPair> &post_c,
            const std::vector<int> &outside_post_vars,
            const int original_op_id);

    void create_operator(int op_id);

    void create_operators_context_split(int op_id);

    void create_operators_context_split_decoupled(int op_id);

    void create_operators();

    void release_memory();

    void dump() const;
};
}

#endif

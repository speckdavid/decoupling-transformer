#ifndef TASKS_DECOUPLED_ROOT_TASK_H
#define TASKS_DECOUPLED_ROOT_TASK_H

#include "root_task.h"

#include "../task_proxy.h"

#include <unordered_map>

namespace plugins {
class Options;
}

namespace decoupling {
class Factoring;
}

namespace tasks {
/*
  Task transformation that decoupled the search space using derived variables and axioms
*/
class DecoupledRootTask : public RootTask {
    std::shared_ptr<RootTask> original_root_task;
    std::shared_ptr<decoupling::Factoring> factoring;

    bool share_effects;
    std::vector<std::vector<ExplicitEffect>> leaf_to_shareable_effects;
    std::vector<std::vector<ExplicitEffect*>> operator_shared_effects;

    std::unordered_map<int, int> center_var_to_pvar;
    std::unordered_map<int, int> leaf_to_goal_svar;
    std::unordered_map<int, std::unordered_map<int, int>> leaf_lstate_to_pvar;
    std::unordered_map<int, std::unordered_map<int, int>> leaf_lstate_to_svar;

    std::unordered_map<int, std::unordered_map<int, int>> leaf_op_to_svar;

public:
    DecoupledRootTask(const plugins::Options &options);
    virtual ~DecoupledRootTask() override = default;

    void reconstruct_plan_if_necessary(std::vector<OperatorID> &path,
                                       std::vector<State> &states) const override;
    
    virtual int get_num_operator_effects(
        int op_index, bool is_axiom) const override;
    virtual int get_num_operator_effect_conditions(
        int op_index, int eff_index, bool is_axiom) const override;
    virtual FactPair get_operator_effect_condition(
        int op_index, int eff_index, int cond_index, bool is_axiom) const override;
    virtual FactPair get_operator_effect(
        int op_index, int eff_index, bool is_axiom) const override;

    virtual TaskProxy get_task_proxy_for_plan_saving() const {
        // If we run decoupled search, we need the original task to save the reconstructed plan.
        return TaskProxy(*original_root_task);
    }

protected:
    virtual const ExplicitEffect &get_effect(int op_id, int effect_id, bool is_axiom) const override;

    void print_statistics() const;
    void write_sas_file(const std::string file_name) const;

    bool are_initial_states_consistent() const;

    void create_variables();
    void create_mutexes();
    void create_initial_state();
    void create_goal();

    // operators
    void precompute_shareable_effects();
    void set_precondition_of_operator(int op_id, ExplicitOperator &op);
    void set_effect_of_operator(int op_id, ExplicitOperator &op);
    void create_operator(int op_id);
    void create_operators();

    // axioms
    void create_frame_axioms();
    void create_goal_axioms();
    void create_precondition_axioms();
    void create_leaf_only_operator_axioms();
    void create_axioms();
};
}

#endif

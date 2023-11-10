#ifndef TASKS_DECOUPLED_ROOT_TASK_H
#define TASKS_DECOUPLED_ROOT_TASK_H

#include "root_task.h"

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

    std::unordered_map<int, int> center_var_to_pvar;
    std::unordered_map<int, int> leaf_to_goal_svar;
    std::unordered_map<int, std::unordered_map<int, int>> leaf_lstate_to_pvar;
    std::unordered_map<int, std::unordered_map<int, int>> leaf_lstate_to_svar;

    std::unordered_map<int, std::unordered_map<int, int>> leaf_op_to_svar;
    std::unordered_map<int, std::unordered_map<int, std::vector<FactPair>>> leaf_op_to_pre;

public:
    DecoupledRootTask(const plugins::Options &options);
    virtual ~DecoupledRootTask() override = default;

protected:
    void print_statistics() const;
    void write_sas_file(const std::string file_name) const;

    bool are_initial_states_consistent() const;

    void create_variables();
    void create_mutexes();
    void create_initial_state();
    void create_goal();

    // operators
    void set_precondition_of_operator(int op_id, ExplicitOperator& op);
    void set_effect_of_operator(int op_id, ExplicitOperator& op);
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

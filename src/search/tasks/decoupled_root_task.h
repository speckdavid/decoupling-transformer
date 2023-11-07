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
    std::unordered_map<int, std::unordered_map<int, int>> leaf_lstate_to_pvar;
    std::unordered_map<int, std::unordered_map<int, int>> leaf_lstate_to_dvar;

public:
    DecoupledRootTask(const plugins::Options &options);
    virtual ~DecoupledRootTask() override = default;

protected:
    void create_variables();
    void create_mutexes();
    void create_initial_state();
    void create_goal();
    void create_operators();
    void create_axioms();
};
}

#endif

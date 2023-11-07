#ifndef TASKS_DECOUPLED_ROOT_TASK_H
#define TASKS_DECOUPLED_ROOT_TASK_H

#include "root_task.h"

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
    std::shared_ptr<decoupling::Factoring> factoring;

public:
    DecoupledRootTask(const plugins::Options &options, const std::shared_ptr<AbstractTask> &abstract_task);
    virtual ~DecoupledRootTask() override = default;
};
}

#endif

#ifndef TASKS_DECOUPLED_TASK_H
#define TASKS_DECOUPLED_TASK_H

#include "delegating_task.h"

namespace plugins {
class Options;
}

namespace extra_tasks {
/*
  Task transformation that decoupled the search space using derived variables and axioms
*/
class DecoupledTask : public tasks::DelegatingTask {

public:
    DecoupledTask(const std::shared_ptr<AbstractTask> &parent);
    virtual ~DecoupledTask() override = default;

};
}

#endif

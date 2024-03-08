#ifndef TASKS_DECOUPLED_PLAN_RECONSTRUCTION_TASK_H
#define TASKS_DECOUPLED_PLAN_RECONSTRUCTION_TASK_H

#include "root_task.h"

#include "../task_proxy.h"

namespace plugins {
class Options;
}

namespace decoupling {
class Factoring;
}

namespace tasks {
/*
  Task transformation that reconstructs a plan (sas_plan) for a given problem (output.sas) and factoring (factoring.txt)
*/
class DecoupledPlanReconstructionTask : public RootTask {
protected:
    std::shared_ptr<decoupling::Factoring> factoring;

public:
    DecoupledPlanReconstructionTask(const plugins::Options &options);
    virtual ~DecoupledPlanReconstructionTask() override = default;
};
}

#endif

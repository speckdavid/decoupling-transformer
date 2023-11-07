#include "decoupled_root_task.h"

#include "../decoupling/factoring.h"

#include "../plugins/plugin.h"
#include "../task_utils/task_dump.h"

#include <iostream>
#include <fstream>

using namespace std;

namespace tasks {
DecoupledRootTask::DecoupledRootTask(const plugins::Options &options)
    : RootTask(),
      original_root_task(dynamic_pointer_cast<RootTask>(tasks::g_root_task)),
      factoring(options.get<shared_ptr<decoupling::Factoring>>("factoring")) {
    // Is this a deep copy?
    variables = original_root_task->variables;
    // mutexes = root_task->mutexes;
    operators = original_root_task->operators;
    axioms = original_root_task->axioms;
    initial_state_values = original_root_task->initial_state_values;
    goals = original_root_task->goals;

    factoring->compute_factoring();
}

class DecoupledRootTaskFeature : public plugins::TypedFeature<AbstractTask, DecoupledRootTask> {
public:
    DecoupledRootTaskFeature() : TypedFeature("decoupled") {
        document_title("Decoupled task");
        document_synopsis(
            "A decoupled transformation of the root task.");

        add_option<shared_ptr<decoupling::Factoring>>("factoring",
                                                      "method that computes the factoring");
    }

    virtual shared_ptr<DecoupledRootTask> create_component(const plugins::Options &options, const utils::Context &) const override {
        return make_shared<DecoupledRootTask>(options);
    }
};

static plugins::FeaturePlugin<DecoupledRootTaskFeature> _plugin;
}

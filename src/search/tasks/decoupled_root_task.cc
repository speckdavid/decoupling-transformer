#include "decoupled_root_task.h"

#include "../decoupling/factoring.h"

#include "../plugins/plugin.h"
#include "../task_utils/task_dump.h"

#include <iostream>
#include <fstream>

using namespace std;

namespace tasks {
DecoupledRootTask::DecoupledRootTask(const plugins::Options &options,
                                     const shared_ptr<AbstractTask> &abstract_task)
    : RootTask(),
      factoring(options.get<shared_ptr<decoupling::Factoring>>("factoring")) {
    shared_ptr<RootTask> root_task = dynamic_pointer_cast<RootTask>(abstract_task);
    if (!root_task) {
        cerr << "Expected a root task as input to decoupled root task." << endl;
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    }
    
    // Is this a deep copy?
    variables = root_task->variables;
    // mutexes = root_task->mutexes;
    operators = root_task->operators;
    axioms = root_task->axioms;
    initial_state_values = root_task->initial_state_values;
    goals = root_task->goals;

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
        return make_shared<DecoupledRootTask>(options, g_root_task);
    }
};

static plugins::FeaturePlugin<DecoupledRootTaskFeature> _plugin;
}

#include "decoupled_task.h"

#include "factoring.h"

using namespace std;

namespace extra_tasks {
DecoupledTask::DecoupledTask(const shared_ptr<AbstractTask> &parent)
    : DelegatingTask(parent) {
        Factoring factoring(parent);
}
}

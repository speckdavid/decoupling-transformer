#include "factoring.h"
#include "../utils/logging.h"

using namespace std;

namespace extra_tasks {
Factoring::Factoring(const std::shared_ptr<AbstractTask> &task) : task(task) {
    utils::g_log << "WARNING: THIS DECOUPLING ONLY WORKS FOR A SPECIALIZED TRANSPORT DOMAIN" << endl;
    create_factoring();
}

void Factoring::create_factoring() {
    for (int var = 0; var < task->get_num_variables(); ++var) {
        utils::g_log << task->get_variable_name(var) << endl;
    }
}
}

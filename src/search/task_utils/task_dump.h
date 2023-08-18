#ifndef TASK_UTILS_TASK_DUMP_H
#define TASK_UTILS_TASK_DUMP_H

#include "../abstract_task.h"

#include <ostream>

namespace task_dump {
extern void dump_as_SAS(const AbstractTask &task, std::ostream &os);
}

#endif

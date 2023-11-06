#ifndef TASKS_FACTORING

#include "../abstract_task.h"

#include <memory>

namespace extra_tasks {
/*
  Dummy Factoring Class to prototype decoupled search
*/

class Factoring {
    std::shared_ptr<AbstractTask> task;
    std::vector<std::vector<int>> leafs;

    void create_factoring();
public:
    Factoring(const std::shared_ptr<AbstractTask> &task);

    int get_leaf_id(int var_id) const;
};
}

#endif

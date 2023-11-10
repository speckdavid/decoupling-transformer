#ifndef DECOUPLING_LEAF_STATE
#define DECOUPLING_LEAF_STATE

#include "leaf_state_id.h"

#include <vector>

// For documentation on classes relevant to storing and working with registered
// states see the file leaf_state_registry.h.

class AbstractTask;
class VariableProxy;

namespace decoupling {

class Factoring;

class LeafState {
    friend class LeafStateRegistry;

    const AbstractTask *task;
    const Factoring *factoring;

    LeafStateID id;

    // Values for vars are maintained in a packed state and accessed on demand.
    std::vector<int> values;

    // Only used by the leaf state registry.
    LeafState(const AbstractTask &task,
              const Factoring &factoring,
              const std::vector<int> &values,
              LeafStateID id);

public:

    LeafStateID get_id() const {
        return id;
    }

    int operator[](size_t index) const;

    int operator[](VariableProxy var) const;

    std::string get_pddl() const;
};
}
#endif

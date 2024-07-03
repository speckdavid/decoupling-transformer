#ifndef DECOUPLING_LEAF_STATE_H
#define DECOUPLING_LEAF_STATE_H

#include "leaf_state_id.h"

#include <vector>

// For documentation on classes relevant to storing and working with registered
// states see the file leaf_state_registry.h.

class AbstractTask;
class VariableProxy;

namespace decoupling {

class Factoring;
class LeafStateRegistry;

class LeafState {
    friend class LeafStateRegistry;

    const AbstractTask *task;
    const Factoring *factoring;
    const LeafStateRegistry *registry;

    const LeafStateID id;

        // Only used by the leaf state registry.
    LeafState(const AbstractTask &task,
              const Factoring &factoring,
              const LeafStateRegistry *registry,
              LeafStateID id);

public:

    LeafStateID get_id() const {
        return id;
    }

    int operator[](int var) const;

    int operator[](VariableProxy var) const;

    std::string get_pddl() const;
};
}
#endif

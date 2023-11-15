#include "leaf_state.h"

#include "factoring.h"
#include "leaf_state_registry.h"

#include <cassert>
#include <iostream>

using namespace std;

namespace decoupling {
LeafState::LeafState(const AbstractTask &task,
                     const Factoring &factoring,
                     const LeafStateRegistry *registry,
                     LeafStateID id)
    : task(&task),
      factoring(&factoring),
      registry(registry),
      id(id) {
    assert(id != LeafStateID::no_state);
    assert(id.hash() != LeafStateHash::MAX);
    assert(id.get_factor() != FactorID::CENTER);
}

int LeafState::operator[](int var) const {
    assert(factoring->get_factor(var) == id.get_factor());
    return registry->leaf_states[id.get_factor()][id.hash()][factoring->get_id_in_factor(var)];
}

int LeafState::operator[](VariableProxy var) const {
    return (*this)[var.get_id()];
}

string LeafState::get_pddl() const {
    string name = "LeafState, id=" + to_string(id.hash()) + ", leaf=" + to_string(id.get_factor());
    for (int var : factoring->get_leaf(id.get_factor())) {
        name += "; ";
        name += task->get_fact_name(FactPair(var, (*this)[var]));
    }
    return name;
}
}
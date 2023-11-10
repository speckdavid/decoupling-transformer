#include "leaf_state.h"

#include "factoring.h"

#include <cassert>
#include <iostream>

using namespace std;

namespace decoupling {
LeafState::LeafState(const AbstractTask &task,
                     const Factoring &factoring,
                     const vector<int> &values,
                     LeafStateID id)
    : task(&task),
      factoring(&factoring),
      id(id),
      values(values){
    assert(id != LeafStateID::no_state);
    assert(id.hash() != LeafStateHash::MAX);
    assert(id.get_factor() != FactorID::CENTER);
}

int LeafState::operator[](size_t index) const {
    assert(factoring->get_factor(index) == id.get_factor());
    return values[factoring->get_id_in_factor(index)];
}

int LeafState::operator[](VariableProxy var) const {
    assert(factoring->get_factor(var.get_id()) == id.get_factor());
    return values[factoring->get_id_in_factor(var.get_id())];
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
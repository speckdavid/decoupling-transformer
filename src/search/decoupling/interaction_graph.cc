#include "interaction_graph.h"

#include <algorithm>
#include <cassert>

using namespace std;

namespace decoupling {
void InteractionGraph::add_dependency(FactorID from, FactorID to) {
    assert(from != to);
    size_t from_id = from;
    if (from == FactorID::CENTER) {
        from_id = num_leaves;
        ifork = false;
    }
    size_t to_id = to;
    if (to == FactorID::CENTER) {
        to_id = num_leaves;
        fork = false;
    }
    if (from != FactorID::CENTER && to != FactorID::CENTER) {
        fork = false;
        ifork = false;
        strict_star = false;
    }
    if (find(successors[from_id].begin(), successors[from_id].end(), to) == successors[from_id].end()) {
        successors[from_id].push_back(to);
    }
    if (find(predecessors[to_id].begin(), predecessors[to_id].end(), from) == predecessors[to_id].end()) {
        predecessors[to_id].push_back(from);
    }
}
}

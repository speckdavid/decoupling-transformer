#ifndef DECOUPLING_INTERACTION_GRAPH
#define DECOUPLING_INTERACTION_GRAPH

#include "leaf_state_id.h"

#include <vector>

namespace decoupling {
class InteractionGraph {
    friend class Factoring;

    // center factor is last

private:
    bool fork, ifork, strict_star;
    int num_leaves;
    std::vector<std::vector<FactorID> > successors;
    std::vector<std::vector<FactorID> > predecessors;

    void add_dependency(FactorID from, FactorID to);

public:
    explicit InteractionGraph(int num_leaves):
            fork(true), ifork(true), strict_star(true), num_leaves(num_leaves) {
        successors.resize(num_leaves+1);
        predecessors.resize(num_leaves+1);
    }

    const std::vector<FactorID> &get_successors(FactorID factor) const {
        if (factor == FactorID::CENTER){
            return successors.back();
        }
        return successors[factor];
    }

    const std::vector<FactorID> &get_predecessors(FactorID factor) const {
        if (factor == FactorID::CENTER){
            return predecessors.back();
        }
        return predecessors[factor];
    }

    bool is_fork_leaf(FactorID leaf) const {
        return successors[leaf].empty();
    }

    bool is_ifork_leaf(FactorID leaf) const {
        return predecessors[leaf].empty();
    }

    bool is_fork() const {
        return fork;
    }

    bool is_ifork() const {
        return ifork;
    }

    bool is_strict_star() const {
        return strict_star;
    }
};
}
#endif

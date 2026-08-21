#ifndef DECOUPLING_INTERACTION_GRAPH_H
#define DECOUPLING_INTERACTION_GRAPH_H

#include "leaf_state_id.h"

#include <cassert>
#include <vector>

namespace decoupling {
class InteractionGraph {
    friend class Factoring;

    // center factor is last

private:
    bool fork, ifork, strict_star;
    int num_leaves;
    std::vector<std::vector<FactorID>> successors;
    std::vector<std::vector<FactorID>> predecessors;

    void add_dependency(FactorID from, FactorID to);

public:
    explicit InteractionGraph(int num_leaves) :
        fork(true), ifork(true), strict_star(true), num_leaves(num_leaves) {
        successors.resize(num_leaves + 1);
        predecessors.resize(num_leaves + 1);
    }

    const std::vector<FactorID> &get_successors(FactorID factor) const {
        if (factor == FactorID::CENTER) {
            return successors.back();
        }
        return successors[factor];
    }

    const std::vector<FactorID> &get_predecessors(FactorID factor) const {
        if (factor == FactorID::CENTER) {
            return predecessors.back();
        }
        return predecessors[factor];
    }

    /*
      A fork leaf must not influence any other factor, i.e. no operator of the
      center or of another leaf may have a precondition on it, and it must not
      be influenced by another leaf, i.e. all operators affecting it may only
      have preconditions on the leaf itself and on the center.
    */
    bool is_fork_leaf(FactorID leaf) const {
        assert(leaf != FactorID::CENTER);
        if (!successors[leaf].empty()) {
            return false;
        }
        if (predecessors[leaf].size() > 1 || 
                (predecessors[leaf].size() == 1 && predecessors[leaf][0] != FactorID::CENTER)) {
            return false;
        }
        return true;
    }

    /*
      Symmetrically, an inverted-fork leaf must not be influenced by any other
      factor, i.e. all operators affecting it may only have preconditions on the
      leaf itself, and it must not influence another leaf, i.e. all operators
      with a precondition on it may only affect the leaf itself and the center.
    */
    bool is_ifork_leaf(FactorID leaf) const {
        assert(leaf != FactorID::CENTER);
        if (!predecessors[leaf].empty()) {
            return false;
        }
        if (successors[leaf].size() > 1 ||
                (successors[leaf].size() == 1 && successors[leaf][0] != FactorID::CENTER)) {
            return false;
        }
        return true;
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

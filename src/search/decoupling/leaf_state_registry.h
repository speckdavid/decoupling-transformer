#ifndef DECOUPLING_LEAF_STATE_REGISTRY_H
#define DECOUPLING_LEAF_STATE_REGISTRY_H

#include "leaf_state.h"
#include "leaf_state_id.h"
#include "../task_proxy.h"

#include "../utils/hash.h"

#include <unordered_map>
#include <vector>
#include <memory>


namespace decoupling {

class LeafStateRegistry {
    friend class LeafState;

    std::shared_ptr<AbstractTask> task;
    std::shared_ptr<Factoring> factoring;

    std::vector<utils::HashMap<std::vector<int>, LeafStateHash>> registered_leaf_states;
    std::vector<std::vector<std::vector<int>>> leaf_states;

    LeafStateHash insert_id_or_pop_leaf_state(FactorID factor);

public:
    LeafStateRegistry(const std::shared_ptr<AbstractTask> &task,
                      const std::shared_ptr<Factoring> &factoring);

    LeafState get_leaf_state(LeafStateHash id, FactorID factor) const;

    LeafStateHash get_leaf_state_hash(const std::vector<int> &facts, FactorID factor);

    LeafStateHash get_successor_leaf_state_hash(const LeafState &predecessor, OperatorProxy op);

    size_t size(FactorID factor) const {
        return registered_leaf_states[factor].size();
    }
};
}
#endif

#include "leaf_state_registry.h"

#include "factoring.h"

using namespace std;

namespace decoupling {
LeafStateRegistry::LeafStateRegistry(const std::shared_ptr<AbstractTask> &task,
                                     const std::shared_ptr<Factoring> &factoring) :
            task(task),
            factoring(factoring) {
    registered_leaf_states.resize(factoring->get_num_leaves());
    leaf_states.resize(factoring->get_num_leaves());

    // add initial leaf states
    vector<int> init_state_data = task->get_initial_state_values();
    for (FactorID leaf(0); leaf < registered_leaf_states.size(); ++leaf){
        vector<int> lstate(factoring->get_leaf(leaf).size(), -1);
        for (int var : factoring->get_leaf(leaf)){
            lstate[factoring->get_id_in_factor(var)] = init_state_data[var];
        }
        registered_leaf_states[leaf].emplace(lstate, 0);
        leaf_states[leaf].push_back(lstate);
    }
}

LeafStateHash LeafStateRegistry::insert_id_or_pop_leaf_state(FactorID factor) {
    const auto [iterator, is_new] = registered_leaf_states[factor].try_emplace(leaf_states[factor].back(),
                                                                               leaf_states[factor].size() - 1);
    if (!is_new) {
        // is known leaf state
        leaf_states[factor].pop_back();
    }
    return iterator->second;
}

LeafState LeafStateRegistry::get_leaf_state(LeafStateHash id, FactorID factor) const {
    return LeafState(*task, *factoring, leaf_states[factor][id], LeafStateID(id, factor));
}

LeafStateHash LeafStateRegistry::get_successor_leaf_state_hash(const LeafState &predecessor, OperatorProxy op) {
    FactorID leaf = predecessor.get_id().get_factor();
    assert(leaf != FactorID::CENTER);
    leaf_states[leaf].push_back(predecessor.values);
    if (factoring->has_eff_on_leaf(op.get_id(), leaf)) {
        for (EffectProxy eff: op.get_effects()) {
            assert(eff.get_conditions().empty());
            if (factoring->get_factor(eff.get_fact().get_variable().get_id()) == leaf) {
                int id_in_leaf = factoring->get_id_in_factor(eff.get_fact().get_variable().get_id());
                leaf_states[leaf].back()[id_in_leaf] = eff.get_fact().get_value();
            }
        }
    }
    return insert_id_or_pop_leaf_state(leaf);
}

LeafStateHash LeafStateRegistry::get_leaf_state_hash(const vector<int> &facts, FactorID factor) {
    assert(factoring->get_leaf(factor).size() == facts.size());
    assert(registered_leaf_states[factor].count(facts) > 0);
    return registered_leaf_states[factor].find(facts)->second;
}
}
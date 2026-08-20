#ifndef DECOUPLING_SIMULATION_RELATION
#define DECOUPLING_SIMULATION_RELATION

#include "leaf_state_id.h"
#include "leaf_state_space.h"

#include <cassert>
#include <vector>
#include <string>
#include <iostream>


namespace utils {
class CountdownTimer;
}
class OperatorID;
class PreconditionsProxy;

namespace decoupling {
class SimulationRelation {

private:
    std::shared_ptr<AbstractTask> task;
    TaskProxy task_proxy;
    std::shared_ptr<Factoring> factoring;
    utils::LogProxy &log;
    LeafStateSpace &leaf_state_space;

    /*
      Label dominance in the center. Whether an operator dominates another one
      only depends on the leaf it affects, its cost, and its center
      preconditions. We therefore group the operators of every fork leaf into
      labels by (cost, center preconditions) and only store the dominance
      relation between the labels of each leaf. This is a block-diagonal
      matrix, where every block is usually tiny compared to the number of
      operators of the leaf. Dominance is only ever checked between operators
      of the same fork leaf, so no cross-leaf entries are needed.
    */
    static constexpr int NO_LABEL = -1;

    /*
      For each operator, its label within the leaf it affects. NO_LABEL for
      global operators and operators of non-fork leaves. Global operators do
      occur as labels of fork-leaf transitions, but they never dominate and are
      never dominated, as they can have center effects.
    */
    std::vector<int> op_to_label;

    // For each leaf, its number of labels (0 for non-fork leaves).
    std::vector<int> num_labels;

    // For each leaf, row-major num_labels x num_labels bit matrix; entry
    // (l, l2) is true iff label l is dominated by label l2 in the center.
    std::vector<std::vector<bool>> label_dominated_by;

    std::vector<std::vector<std::vector<bool>>> relation;


    // both fact vectors must be sorted; true iff pre is a subset of pre2
    bool center_precondition_dominance(const std::vector<FactPair> &pre,
                                       const std::vector<FactPair> &pre2) const;

    // true iff op is dominated by op2; both must have an effect on leaf factor
    inline bool op_dominated_by(FactorID factor, OperatorID op, OperatorID op2) const {
        int label = op_to_label[op.get_index()];
        int label2 = op_to_label[op2.get_index()];
        if (label == NO_LABEL || label2 == NO_LABEL){
            // at least one of them is a global operator
            return false;
        }
        assert(label < num_labels[factor] && label2 < num_labels[factor]);
        return label_dominated_by[factor][label * num_labels[factor] + label2];
    }

    void compute_label_dominance();

    void compute_simulation(FactorID factor, const utils::CountdownTimer &timer);

    inline bool simulates(FactorID factor, LeafStateHash s, LeafStateHash t) const {
        return relation[factor][s][t];
    }

    inline bool similar(FactorID factor, LeafStateHash s, LeafStateHash t) const {
        return relation[factor][s][t] && relation[factor][t][s];
    }

    inline void remove(FactorID factor, LeafStateHash s, LeafStateHash t) {
        relation[factor][s][t] = false;
    }

    inline const std::vector<std::vector<bool>> & get_relation(FactorID factor) const {
        return relation[factor];
    }

    size_t num_equivalences(FactorID factor) const;

    size_t num_simulations(FactorID factor, bool ignore_equivalences) const;

    size_t num_states(FactorID factor) const {
        return relation[factor].size();
    }

    void dump(FactorID factor) const;

    std::string get_name(LeafStateHash id, FactorID factor) const;

public:

    SimulationRelation(const std::shared_ptr<AbstractTask> &task,
                       const std::shared_ptr<Factoring> &factoring,
                       utils::LogProxy &log,
                       LeafStateSpace &leaf_state_space,
                       int time_limit);

    void statistics() const;

    void perform_leaf_irrelevance_pruning(bool prune_bwd_graph, bool only_remove_states);

};
}
#endif
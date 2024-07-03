#ifndef DECOUPLING_LEAF_STATE_SPACE_H
#define DECOUPLING_LEAF_STATE_SPACE_H

#include "leaf_state_id.h"

#include "../operator_id.h"
#include "../task_proxy.h"

#include "../utils/logging.h"

#include <memory>
#include <vector>


namespace decoupling {
class Factoring;
class LeafState;
class LeafStateRegistry;

class LeafStateSpace {
    friend class PathPrices;
    friend class SimulationRelation;

    utils::LogProxy &log;

    std::shared_ptr<Factoring> factoring;
    std::shared_ptr<AbstractTask> task;
    std::shared_ptr<LeafStateRegistry> state_registry;

    TaskProxy task_proxy;

    int num_leaves;

    std::vector<std::vector<bool>> is_a_leaf_goal_state;

    std::vector<std::vector<LeafStateHash>> leaf_goal_states;

    std::vector<std::vector<std::vector<std::pair<OperatorID, LeafStateHash>>>> leaf_state_successors;

    std::vector<std::vector<std::vector<std::pair<OperatorID, LeafStateHash>>>> leaf_state_predecessors;

    std::vector<bool> is_leaf_state_space_scc;


    bool has_center_precondition(OperatorProxy op) const;

    void get_applicable_ops(const LeafState &lstate, std::vector<OperatorID> &applicable_ops) const;

    bool check_is_goal_state(const LeafState &lstate);

    void check_leaf_invertibility(FactorID leaf,
                                  const std::vector<std::vector<int>> &leaf_only_state_space_graph);

    void build_leaf_state_space(FactorID leaf,
                                bool compute_leaf_invertibility);

    /*
      this builds all leaf state spaces and stores them
    */
    void build_leaf_state_spaces(bool compute_leaf_invertibility,
                                 bool prune_fork_leaf_state_spaces);

public:
    LeafStateSpace(const std::shared_ptr<Factoring> &factoring,
                   const std::shared_ptr<AbstractTask> &task,
                   utils::LogProxy &log,
                   bool compute_leaf_invertibility,
                   bool prune_fork_leaf_state_spaces);

    LeafState get_leaf_state(LeafStateHash id, FactorID leaf) const;

    bool is_leaf_goal_state(LeafStateHash id, FactorID factor) const {
        assert(is_a_leaf_goal_state.size() > factor && is_a_leaf_goal_state[factor].size() > id);
        return is_a_leaf_goal_state[factor][id];
    }

    const std::vector<LeafStateHash> &get_leaf_goal_states(FactorID factor) const {
        assert(leaf_goal_states.size() > factor);
        return leaf_goal_states[factor];
    }

    bool is_applicable(LeafStateHash id, FactorID leaf, OperatorProxy op) const;

    bool is_applicable(LeafStateHash id, FactorID leaf, OperatorID op_id) const;

    std::vector<int> get_predecessors(FactorID leaf, LeafStateHash lstate, OperatorID op_id) const {
        std::vector<int> predecessors;
        for (auto [p_op_id, predecessor] : leaf_state_predecessors[leaf][lstate]) {
            if (op_id == p_op_id) {
                predecessors.push_back(predecessor);
            }
        }
        return predecessors;
    }

    size_t get_num_states(FactorID leaf) const {
        assert(leaf != FactorID::CENTER);
        return is_a_leaf_goal_state[leaf].size();
    }

    bool is_leaf_invertible(FactorID leaf) const {
        assert(is_leaf_state_space_scc.size() > leaf);
        return is_leaf_state_space_scc[leaf];
    }

    std::string get_name(LeafStateHash id, FactorID leaf) const;
};
}
#endif

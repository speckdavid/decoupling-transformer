#ifndef DECOUPLING_PLAN_RECONSTRUCTION_H
#define DECOUPLING_PLAN_RECONSTRUCTION_H

#include "leaf_state_id.h"

#include "../operator_id.h"

#include <vector>

class AbstractTask;
class OperatorProxy;
class State;
class TaskProxy;

namespace decoupling {
class Factoring;
class LeafStateSpace;

struct PathPriceInfo {
    int price;
    OperatorID generating_op;
    LeafStateHash predecessor;
    bool is_new;

    PathPriceInfo() : price(-1),
                      generating_op(OperatorID::no_operator),
                      predecessor(LeafStateHash::MAX),
                      is_new(false) {}

    PathPriceInfo(int price, OperatorID op, LeafStateHash predecessor, bool is_new = false)
        : price(price), generating_op(op), predecessor(predecessor), is_new(is_new) {}

    void reset_generating_op() {
        is_new = false;
        generating_op = OperatorID::no_operator;
    }

    void dump(const AbstractTask &task) const;
};


class PathPrices {
    std::vector<int> number_states;

    std::vector<int> goal_costs;

    // path information
    std::vector<std::vector<PathPriceInfo>> paths;

    // return true, if the LeafState was actually added
    bool add_state(LeafStateHash id, FactorID leaf, int cost,
                   OperatorID generating_op, LeafStateHash predecessor,
                   bool is_goal_state);

    explicit PathPrices(int num_leaves);

    PathPrices(const PathPrices &cpg);

    int get_number_states(FactorID leaf) const {
        return number_states[leaf];
    }

    bool has_leaf_state(LeafStateHash id, FactorID leaf) const {
        return id < paths[leaf].size() && paths[leaf][id].price != -1;
    }

    int get_cost_of_state(LeafStateHash id, FactorID leaf) const {
        return paths[leaf][id].price;
    }

    const PathPriceInfo &get_path_info(LeafStateHash id, FactorID factor) const;

    void update(const State &new_center_state,
                const TaskProxy &task_proxy,
                const Factoring &factoring,
                const LeafStateSpace &leaf_state_space);

    void apply_global_op_to_leaves(const PathPrices &old_cpg,
                                   OperatorProxy op,
                                   const Factoring &factoring,
                                   const LeafStateSpace &leaf_state_space);

public:

    static void insert_leaf_actions(const AbstractTask &task,
                                    const Factoring &factoring,
                                    const LeafStateSpace &leaf_state_space,
                                    std::vector<OperatorID> &ops,
                                    std::vector<State> &states);
};
}
#endif

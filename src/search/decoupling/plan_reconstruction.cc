#include "plan_reconstruction.h"

#include "factoring.h"
#include "leaf_state_registry.h"
#include "../utils/system.h"

#include <limits>
#include <iostream>
#include <numeric>

using namespace std;

namespace decoupling {
void PathPriceInfo::dump(const AbstractTask &task) const {
    cout << "PathPriceTagInfo" << endl;
    if (generating_op != OperatorID::no_operator){
        cout << "op: " << task.get_operator_name(generating_op.get_index(), false) << endl;
    } else {
        cout << "initial state fact" << endl;
    }
    if (is_new){
        cout << "is new" << endl;
    }
}

PathPrices::PathPrices(int num_leaves) {
    number_states.resize(num_leaves, 0);
    goal_costs.resize(num_leaves, -1);
    paths.resize(num_leaves);
}

PathPrices::PathPrices(const PathPrices &other) {
    number_states = other.number_states;
    goal_costs = other.goal_costs;
    paths = other.paths;
    for (FactorID factor(0); factor < paths.size(); ++factor){
        for (LeafStateHash state(0); state < paths[factor].size(); ++state){
            paths[factor][state].reset_generating_op();
        }
    }
}

bool PathPrices::add_state(LeafStateHash id, FactorID leaf, int cost,
                           OperatorID generating_op, LeafStateHash predecessor,
                           bool is_goal_state) {
    bool added = false;
    if (id >= paths[leaf].size()){
        added = true;
        paths[leaf].resize(id + 1);
        number_states[leaf]++;
    } else if (paths[leaf][id].price == -1) {
        added = true;
        number_states[leaf]++;
    } else if (paths[leaf][id].price > cost){
        added = true;
    }
    if (added){
        paths[leaf][id].price = cost;
        paths[leaf][id].generating_op = generating_op;
        paths[leaf][id].predecessor = predecessor;
        paths[leaf][id].is_new = true;
        if (is_goal_state && (goal_costs[leaf] == -1 || goal_costs[leaf] > cost)){
            goal_costs[leaf] = cost;
        }
    }
    return added;
}

const PathPriceInfo& PathPrices::get_path_info(LeafStateHash id, FactorID factor) const {
    return paths[factor][id];
}

void PathPrices::update(const State &base_state,
                        const TaskProxy &task_proxy,
                        const Factoring &factoring,
                        const LeafStateSpace &leaf_state_space) {

    for (FactorID leaf(0); leaf < factoring.get_num_leaves(); ++leaf){

        if (factoring.is_fork_leaf(leaf) && !factoring.has_leaf_goal(leaf)){
            // skip fork leaves that don't have a goal or whose goal is already achieved (only satisficing search)
            continue;
        }

        vector<int> best_prices(leaf_state_space.get_num_states(leaf), numeric_limits<int>::max());

        bool change = true;

        while (change) {

            change = false;

            int num_states = 0;
            LeafStateHash id(0);
            while (num_states < number_states[leaf]) {
                if (has_leaf_state(id, leaf)){
                    ++num_states;

                    int cost = get_cost_of_state(id, leaf);
                    if (best_prices[id] <= cost){
                        ++id;
                        continue;
                    }

                    best_prices[id] = cost;

                    for (const auto &[op_id, succ_id] : leaf_state_space.leaf_state_successors[leaf][id]){
                        if (factoring.is_global_operator(op_id.get_index())){
                            continue;
                        }
                        OperatorProxy op = task_proxy.get_operators()[op_id];
                        if (factoring.is_ifork_leaf(leaf) ||
                            factoring.is_center_applicable(base_state, op)){
                            change |= add_state(succ_id, leaf,
                                                cost + task_proxy.get_operators()[op.get_id()].get_cost(),
                                                OperatorID(op.get_id()),
                                                id,
                                                leaf_state_space.is_leaf_goal_state(succ_id, leaf));
                        }
                    }
                }
                ++id;
            }
        }
    }
}

void PathPrices::apply_global_op_to_leaves(const PathPrices &old_cpg,
                                           OperatorProxy op,
                                           const Factoring &factoring,
                                           const LeafStateSpace &leaf_state_space) {
    for (FactorID leaf(0); leaf < factoring.get_num_leaves(); ++leaf){

        bool has_pre_on_factor = factoring.has_pre_on_leaf(OperatorID(op.get_id()), leaf);
        size_t num_effects = factoring.get_num_effects_on_leaf(op, leaf);

        if (!has_pre_on_factor && num_effects == 0){
            // op does not affect leaf at all => copy everything
            number_states[leaf] = old_cpg.number_states[leaf];
            goal_costs[leaf] = old_cpg.goal_costs[leaf];
            paths[leaf] = old_cpg.paths[leaf];
            for (size_t id = 0; id < paths[leaf].size(); ++id){
                paths[leaf][id].reset_generating_op();
            }
            continue;
        }

        if (num_effects == factoring.get_leaf(leaf).size()){
            // center op completely overwrites leaf
            LeafState predecessor = leaf_state_space.get_leaf_state(LeafStateHash(0), leaf);
            LeafStateHash succ_state = leaf_state_space.state_registry->get_successor_leaf_state_hash(predecessor, op);
            add_state(succ_state,
                      leaf,
                      0,
                      OperatorID::no_operator,
                      LeafStateHash::MAX,
                      leaf_state_space.is_leaf_goal_state(succ_state, leaf));
            paths[leaf][succ_state].reset_generating_op();
            continue;
        }

        // check which leaf states satisfy center precondition and apply effects
        size_t number_states = old_cpg.get_number_states(leaf);
        LeafStateHash id(0);
        while (number_states > 0) {
            if (old_cpg.has_leaf_state(id, leaf)){
                --number_states;

                if (!has_pre_on_factor || leaf_state_space.is_applicable(id, leaf, op)){
                    if (num_effects == 0){
                        add_state(id, leaf,
                                  old_cpg.get_cost_of_state(id, leaf),
                                  OperatorID::no_operator,
                                  id,
                                  leaf_state_space.is_leaf_goal_state(id, leaf));
                        paths[leaf][id].reset_generating_op();
                    } else {
                        LeafState state = leaf_state_space.get_leaf_state(id, leaf);
                        LeafStateHash succ_id = leaf_state_space.state_registry->get_successor_leaf_state_hash(state, op);

                        add_state(succ_id, leaf,
                                  old_cpg.get_cost_of_state(id, leaf),
                                  OperatorID::no_operator,
                                  id,
                                  leaf_state_space.is_leaf_goal_state(succ_id, leaf));
                        paths[leaf][succ_id].reset_generating_op();
                    }
                }
            }
            ++id;
        }
    }
}

void PathPrices::insert_leaf_actions(const AbstractTask &task,
                                     const Factoring &factoring,
                                     const LeafStateSpace &leaf_state_space,
                                     vector<OperatorID> &path,
                                     vector<State> &states) {

    TaskProxy task_proxy(task);

    int num_leaves = factoring.get_num_leaves();

    reverse(path.begin(), path.end());
    reverse(states.begin(), states.end());

    // reconstruct leaf paths
    vector<unique_ptr<PathPrices>> price_tags(states.size());

    price_tags[0].reset(new PathPrices(num_leaves));

    for (FactorID leaf(0); leaf < num_leaves; ++leaf){
        price_tags[0]->add_state(LeafStateHash(0),
                                 leaf,
                                 0,
                                 OperatorID::no_operator,
                                 LeafStateHash::MAX,
                                 leaf_state_space.is_leaf_goal_state(LeafStateHash(0), leaf));
    }
    price_tags[0]->update(states[0], task_proxy, factoring, leaf_state_space);

    for (size_t i = 1; i < states.size(); ++i){
        if (!factoring.is_fork_factoring()){
            price_tags[i].reset(new PathPrices(num_leaves));
            price_tags[i]->apply_global_op_to_leaves(*price_tags[i - 1],
                                                     task_proxy.get_operators()[path[i - 1]],
                                                     factoring,
                                                     leaf_state_space);
        } else {
            price_tags[i].reset(new PathPrices(*price_tags[i - 1].get()));
        }
        price_tags[i]->update(states[i], task_proxy, factoring, leaf_state_space);
    }

    // start from goal state
    reverse(path.begin(), path.end());
    reverse(price_tags.begin(), price_tags.end());

    vector<LeafStateHash> marked_leaf_states(num_leaves, LeafStateHash::MAX);

    // mark leaf goal states
    for (FactorID leaf(0); leaf < num_leaves; ++leaf){
        if (factoring.has_leaf_goal(leaf)){
            int min_cost = numeric_limits<int>::max();

            size_t number_states = price_tags[0]->get_number_states(leaf);
            LeafStateHash id(0);
            while (number_states > 0) {
                if (price_tags[0]->has_leaf_state(id, leaf)){
                    --number_states;
                    if (leaf_state_space.is_leaf_goal_state(id, leaf)){
                        int new_cost = price_tags[0]->get_cost_of_state(id, leaf);
                        assert(new_cost >= 0);
                        if (min_cost > new_cost){
                            min_cost = new_cost;
                            marked_leaf_states[leaf] = id;
                        }
                    }
                }
                ++id;
            }
            assert(marked_leaf_states[leaf] != LeafStateHash::MAX);
        }
    }

    vector<OperatorID> decoupled_plan;

    // go through center path and fill in leaf actions
    for (size_t step = 0; step < price_tags.size(); ++step){

        OperatorID op_id(OperatorID::no_operator);
        if (step < path.size()){
            op_id = path[step];
        }

        // iterate over leaf factors
        for (FactorID factor(0); factor < num_leaves; ++factor){

            bool change = true;
            while (change && marked_leaf_states[factor] != LeafStateHash::MAX){
                // backtrack in current CPG to next center action
                change = false;

                const PathPriceInfo &path_info = price_tags[step]->get_path_info(marked_leaf_states[factor], factor);
                if (path_info.is_new){
                    if (path_info.generating_op != OperatorID::no_operator){
                        // leaf action
                        assert(path_info.predecessor != LeafStateHash::MAX);
                        assert(path_info.predecessor != marked_leaf_states[factor]); // not a self-loop
                        assert(factoring.is_leaf_only_operator(path_info.generating_op.get_index()) &&
                               find(factoring.get_leaf_operators(factor).begin(),
                                    factoring.get_leaf_operators(factor).end(),
                                    path_info.generating_op) != factoring.get_leaf_operators(factor).end());

                        marked_leaf_states[factor] = path_info.predecessor;
                        decoupled_plan.push_back(path_info.generating_op);
                        change = true;
                    } else {
                        // marked_leaf_states[factor] is initial leaf state
                        assert(op_id == OperatorID::no_operator);
                        assert(path_info.predecessor == LeafStateHash::MAX);

                        marked_leaf_states[factor] = LeafStateHash::MAX;
                    }
                }
            }

            bool need_compliant_leaf_state = false;

            if (marked_leaf_states[factor] != LeafStateHash::MAX &&
                op_id != OperatorID::no_operator &&
                    factoring.has_eff_on_leaf(op_id, factor)){
                // mimic leaf effects to predecessor decoupled state
                const PathPriceInfo &path_info = price_tags[step]->get_path_info(marked_leaf_states[factor], factor);

                assert(path_info.generating_op == OperatorID::no_operator);

                if (path_info.predecessor != LeafStateHash::MAX){
                    // predecessor leaf state in predecessor decoupled state
                    marked_leaf_states[factor] = path_info.predecessor;
                } else {
                    // center action completely overwrites leaf,
                    // pick any leaf state in predecessor decoupled state
//                    assert(op->get_effects(factor).size() == factoring.get_leaf(factor).size()); TODO

                    need_compliant_leaf_state = true;
                    marked_leaf_states[factor] = LeafStateHash::MAX;
                }
            }

            if (need_compliant_leaf_state || (op_id != OperatorID::no_operator &&
                    factoring.has_eff_on_leaf(op_id, factor))){
                // collect and mark leaf precondition of center op
                assert(marked_leaf_states[factor] == LeafStateHash::MAX ||
                       leaf_state_space.is_applicable(marked_leaf_states[factor], factor, op_id));

                if (marked_leaf_states[factor] == LeafStateHash::MAX){
                    int best_price = numeric_limits<int>::max();

                    size_t number_states = price_tags[step + 1]->get_number_states(factor);
                    LeafStateHash id(0);
                    while (number_states > 0) {
                        if (price_tags[step + 1]->has_leaf_state(id, factor)){
                            --number_states;
                            if (leaf_state_space.is_applicable(id, factor, op_id)){
                                int cost = price_tags[step + 1]->get_cost_of_state(id, factor);
                                if (cost < best_price){
                                    best_price = cost;
                                    marked_leaf_states[factor] = id;
                                }
                            }
                        }
                        ++id;
                    }

                    assert(marked_leaf_states[factor] != LeafStateHash::MAX);
                }
            }
        }

        if (op_id != OperatorID::no_operator){
            decoupled_plan.push_back(op_id);
        }
    }

#ifndef NDEBUG
    for (int i = 0; i < num_leaves; ++i){
        // all leaves are in initial state
        assert(marked_leaf_states[i] == LeafStateHash::MAX);
    }
#endif

    path.swap(decoupled_plan);
}
}
#include "leaf_state_space.h"

#include "factoring.h"
#include "leaf_state.h"
#include "leaf_state_registry.h"

#include "../algorithms/sccs.h"

#include <limits>

using namespace std;

namespace decoupling {

LeafStateSpace::LeafStateSpace(const shared_ptr<Factoring> &factoring,
                               const shared_ptr<AbstractTask> &task,
                               utils::LogProxy &log,
                               bool compute_leaf_invertibility,
                               bool prune_fork_leaf_state_spaces) :
    log(log),
    factoring(factoring),
    task(task),
    state_registry(make_shared<LeafStateRegistry>(LeafStateRegistry(task, factoring))),
    task_proxy(TaskProxy(*task)),
    num_leaves(factoring->get_leaves().size()) {

    leaf_state_successors.resize(num_leaves);
    leaf_state_predecessors.resize(num_leaves);

    is_leaf_state_space_scc.resize(num_leaves);

    is_a_leaf_goal_state.resize(num_leaves);
    leaf_goal_states.resize(num_leaves);

    build_leaf_state_spaces(compute_leaf_invertibility, prune_fork_leaf_state_spaces);
}

LeafState LeafStateSpace::get_leaf_state(LeafStateHash id, FactorID leaf) const {
    return state_registry->get_leaf_state(id, leaf);
}

bool LeafStateSpace::has_center_precondition(OperatorProxy op) const {
    for (FactProxy pre : op.get_preconditions()){
        if (factoring->get_factor(pre.get_variable().get_id()) == FactorID::CENTER){
            return true;
        }
    }
    return false;
}

bool LeafStateSpace::is_applicable(LeafStateHash id, FactorID leaf, OperatorProxy op) const {
    assert(leaf != FactorID::CENTER);
    LeafState lstate = state_registry->get_leaf_state(id, leaf);
    for (FactProxy pre : op.get_preconditions()){
        if (factoring->get_factor(pre.get_variable().get_id()) != leaf){
            continue;
        }
        if (lstate[pre.get_variable().get_id()] != pre.get_value()){
            return false;
        }
    }
    return true;
}

void LeafStateSpace::get_applicable_ops(const LeafState &lstate,
                                        vector<OperatorID> &applicable_ops) const {
    FactorID leaf = lstate.get_id().get_factor();
    assert(leaf != FactorID::CENTER);

    for (OperatorID op_id : factoring->get_leaf_operators(leaf)){
        OperatorProxy op = task_proxy.get_operators()[op_id];
        bool applicable = true;
        for (FactProxy pre : op.get_preconditions()){
            if (lstate[pre.get_variable()] != pre.get_value()){
                applicable = false;
                break;
            }
        }
        if (applicable){
            applicable_ops.emplace_back(op.get_id());
        }
    }
}

bool LeafStateSpace::check_is_goal_state(const LeafState &lstate) {
    FactorID leaf = lstate.get_id().get_factor();
    assert(leaf != FactorID::CENTER);
    for (FactPair leaf_goal : factoring->get_leaf_goals(leaf)){
        assert(leaf == factoring->get_factor(leaf_goal.var));
        if (lstate[leaf_goal.var] != leaf_goal.value){
            return false;
        }
    }
    return true;
}

void LeafStateSpace::check_leaf_invertibility(FactorID leaf,
                                              const vector<vector<int>> &leaf_only_state_space_graph) {
    size_t num_sccs = sccs::compute_maximal_sccs(leaf_only_state_space_graph).size();
    if (num_sccs == 1){
        size_t prod_size = 1;
        for (int var : factoring->get_leaf(leaf)){
            prod_size *= task->get_variable_domain_size(var);
        }
        log << "state space of leaf " << leaf << " is strongly connected via leaf-only actions" << endl;

        is_leaf_state_space_scc[leaf] = true;

        if (prod_size != state_registry->size(leaf)){
            // TODO could do this for all leaves where leaf state space is constructed
            log << "WARNING: not all leaf states for leaf " << leaf << " are reachable"
                 << ", removing non-applicable center actions from successor generator" << endl;
            factoring->remove_never_applicable_global_ops(leaf);
        }
    }
}

void LeafStateSpace::build_leaf_state_space(FactorID factor,
                                            bool compute_leaf_invertibility) {

    vector<bool> closed(1, false); // initial state
    bool change = true;

    vector<vector<int>> leaf_only_state_space_graph;

    while (change) {

        change = false;

        for (LeafStateHash id(0); id < closed.size(); ++id){
            if (id < closed.size()){
                if (closed[id]) {
                    // no need to handle a state twice
                    continue;
                }
            } else {
                closed.resize(id + 1, false);
                if (compute_leaf_invertibility){
                    leaf_only_state_space_graph.resize(id + 1);
                }
            }

            const LeafState curr_leaf_state = state_registry->get_leaf_state(id, factor);
            closed[id] = true;
            check_is_goal_state(curr_leaf_state);

            if (id >= leaf_state_successors[factor].size()) {
                leaf_state_successors[factor].resize(id + 1);
            }
            if (id >= leaf_state_predecessors[factor].size()){
                leaf_state_predecessors[factor].resize(id + 1);
            }

            if (factoring->is_fork_leaf(factor) && is_a_leaf_goal_state[factor][id]){
                // the goal states of fork leaves don't need successors
                continue;
            }

            vector<OperatorID> applicable_ops;
            get_applicable_ops(curr_leaf_state, applicable_ops);

            // apply all applicable_ops to predecessor and store the outcome
            for (OperatorID op_id : applicable_ops){
                OperatorProxy op = task_proxy.get_operators()[op_id];
                bool is_global_op = factoring->is_global_operator(op_id.hash());

                LeafStateHash succ_id = state_registry->get_successor_leaf_state_hash(curr_leaf_state, op);

                if (!is_global_op && id == succ_id){
                    // this would induce a self-loop in the leaf state space
                    // need to keep track of these for global operators
                    continue;
                }

                if (compute_leaf_invertibility &&
                    !has_center_precondition(op) &&
                    !is_global_op){
                    // is leaf action without center precondition
                    leaf_only_state_space_graph[id].push_back(succ_id);
                }

                // this ensures that there is an entry leaf_state_predecessors[factor][lid]
                // even in case succ_id has been generated by a center action
                if (succ_id >= leaf_state_predecessors[factor].size()) {
                    leaf_state_predecessors[factor].resize(succ_id + 1);
                }

                leaf_state_successors[factor][id].emplace_back(op_id, succ_id);
                leaf_state_predecessors[factor][succ_id].emplace_back(op_id, id);
            }
        }
    }

    if (compute_leaf_invertibility){
        check_leaf_invertibility(factor, leaf_only_state_space_graph);
    }
}

void LeafStateSpace::build_leaf_state_spaces(bool compute_leaf_invertibility,
                                             bool prune_fork_leaf_state_spaces) {

    for (FactorID factor(0); factor < num_leaves; ++factor){
        if (factoring->is_fork_leaf(factor) && !factoring->has_leaf_goal(factor)){
            // skip fork-leaf factor without a goal
            continue;
        }
        build_leaf_state_space(factor, compute_leaf_invertibility);
    }

    if (prune_fork_leaf_state_spaces) {
        // TODO integrate this
        log << "ERROR: pruning leaf states spaces is not yet implemented" << endl;
//    pruning->apply_leaf_state_space_pruning();
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }

    size_t min_leaf_factor_size = numeric_limits<size_t>::max();
    double avg_leaf_factor_size = 0;
    size_t max_leaf_factor_size = 0;
    size_t tmp;

    for (FactorID factor(0); factor < num_leaves; ++factor){
        tmp = state_registry->size(factor);
        min_leaf_factor_size = min(min_leaf_factor_size, tmp);
        avg_leaf_factor_size += tmp;
        max_leaf_factor_size = max(max_leaf_factor_size, tmp);
    }

    log << "min reachable leaf factor size "  << min_leaf_factor_size << endl;
    log << "avg reachable leaf factor size "  << (int) (avg_leaf_factor_size/num_leaves) << endl;
    log << "max reachable leaf factor size "  << max_leaf_factor_size << endl;

    for (FactorID factor(0); factor < num_leaves; ++factor){
        if (factoring->has_leaf_goal(factor) && leaf_goal_states[factor].empty()){
            log << "There is a goal that cannot be achieved in factor " << factor << endl;
            for (FactPair goal : factoring->get_leaf_goals(factor)){
                log << task->get_fact_name(goal) << endl;
            }
            exit_with(utils::ExitCode::SEARCH_UNSOLVABLE);
        }
    }
    log << "done building leaf state spaces" << endl;
}

string LeafStateSpace::get_name(LeafStateHash id, FactorID leaf) const {
    return state_registry->get_leaf_state(id, leaf).get_pddl();
}
}

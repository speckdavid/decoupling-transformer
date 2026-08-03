#include "simulation_relation.h"

#include "factoring.h"

#include "../algorithms/priority_queues.h"
#include "../operator_id.h"
#include "../utils/countdown_timer.h"
#include "../utils/timer.h"

using namespace std;

namespace decoupling {
SimulationRelation::SimulationRelation(const std::shared_ptr<AbstractTask> &task,
                                       const std::shared_ptr<Factoring> &factoring,
                                       utils::LogProxy &log,
                                       LeafStateSpace &leaf_state_space,
                                       int time_limit) :
        task(task),
        task_proxy(TaskProxy(*task)),
        factoring(factoring),
        log(log),
        leaf_state_space(leaf_state_space) {

    utils::CountdownTimer timer(time_limit > 0 ? time_limit : numeric_limits<double>::infinity());
    log << "Initializing simulation relation." << endl;

    bool has_fork_leaf = factoring->is_fork_factoring();
    if (!has_fork_leaf){
        for (FactorID leaf(0); leaf < factoring->get_num_leaves(); ++leaf){
            if (factoring->is_fork_leaf(leaf)) {
                has_fork_leaf = true;
                break;
            }
        }
    }
    if (!has_fork_leaf){
        log << "No fork leaf found, skipping irrelevance pruning." << endl;
        return;
    }

    compute_label_dominance();
    if (timer.is_expired()){
        return;
    }
    relation.resize(factoring->get_num_leaves());
    for (FactorID factor(0); factor < relation.size(); ++factor){
        if (factoring->is_fork_leaf(factor)){
            log  << "Computing simulation relation for leaf " << factor << endl;
            compute_simulation(factor, timer);
            if (timer.is_expired()){
                return;
            }
        }
    }
}

void SimulationRelation::perform_leaf_irrelevance_pruning(bool prune_bwd_graph, bool only_remove_states) {
    vector<vector<vector<pair<OperatorID, LeafStateHash> > > > &transition_system_fwd = leaf_state_space.leaf_state_successors;
    vector<vector<vector<pair<OperatorID, LeafStateHash> > > > &transition_system_bwd = leaf_state_space.leaf_state_predecessors;

    vector<vector<vector<pair<OperatorID, LeafStateHash> > > > tmp;
    if (only_remove_states){
        tmp = transition_system_fwd;
    }

    int num_transitions_before = 0, num_transitions_after = 0;
    for (FactorID factor(0); factor < relation.size(); ++factor){
        if (relation[factor].empty()){
            // is a non-fork leaf or the computation timed out
            continue;
        }
        for (LeafStateHash s(0); s  < transition_system_fwd[factor].size(); ++s){
            auto  & trs_s = transition_system_fwd[factor][s];
            num_transitions_before += trs_s.size();
            trs_s.erase(std::remove_copy_if(std::begin(trs_s),
                                            std::end(trs_s), std::begin(trs_s),
                                            [&](pair<OperatorID, LeafStateHash> & tr){
                                                if (simulates(factor, s, tr.second)) return true;
                                                return std::find_if(std::begin(trs_s),
                                                                    std::end(trs_s),
                                                                    [&](pair<OperatorID, LeafStateHash> & tr2){
                                                                        return simulates(factor, tr2.second, tr.second) &&
                                                                            op_dominated_by[tr.first.hash()][tr2.first.hash()]
                                                                            && (!(simulates(factor, tr.second, tr2.second) &&
                                                                                  op_dominated_by[tr2.first.hash()] [tr.first.hash()])
                                                                                || tr.second < tr2.second ||
                                                                                (tr.second == tr2.second &&
                                                                                 tr.first.hash() < tr2.first.hash()));
                                                                    }) != std::end(trs_s);
                                            }),
                        std::end(trs_s));
        }

        // forward reachability analysis; are states reachable after irrelevance pruning
        vector<bool> reachable(transition_system_fwd[factor].size(), false);

        vector<size_t> current(1, 0);

        while (!current.empty()){
            vector<size_t> next;
            for (auto id : current){
                if (reachable[id]){
                    continue;
                }
                reachable[id] = true;
                for (const auto &transition : transition_system_fwd[factor][id]){
                    int t = transition.second;
                    if (!reachable[t]) {
                        next.push_back(t);
                    }
                }
            }
            next.swap(current);
        }

        if (only_remove_states){
            transition_system_fwd[factor].swap(tmp[factor]);
        }

        // remove irrelevant states
        for (LeafStateHash id(0); id < reachable.size(); ++id){
            if (!reachable[id]){
                vector<pair<OperatorID, LeafStateHash> >().swap(transition_system_fwd[factor][id]);
            }
        }

        // (re-)remove transitions entering removed states
        for (LeafStateHash s(0); s  < transition_system_fwd[factor].size(); ++s){
            auto &trs_s = transition_system_fwd[factor][s];

            trs_s.erase(remove_if(std::begin(trs_s), std::end(trs_s),
                    [&](pair<OperatorID, LeafStateHash> & tr){
                return !reachable[tr.second];
            }), std::end(trs_s));

            num_transitions_after += trs_s.size();
        }
    }

    if (prune_bwd_graph) {
        for (FactorID factor(0); factor < relation.size(); ++factor){
            if (relation[factor].empty()){
                // is a non-fork leaf or the computation timed out
                continue;
            }
            for (LeafStateHash s(0); s  < transition_system_bwd[factor].size(); ++s){
                vector<pair<OperatorID, LeafStateHash> >().swap(transition_system_bwd[factor][s]);
            }

            for (LeafStateHash s(0); s  < transition_system_fwd[factor].size(); ++s){
                for (const auto &tr : transition_system_fwd[factor][s]){
                    transition_system_bwd[factor][tr.second].push_back(make_pair(tr.first, s));
                }
            }
        }
    }

    log << "Irrelevance pruning: " << num_transitions_before << " => " << num_transitions_after << " transitions remaining" << endl;
}

bool SimulationRelation::center_precondition_dominance(
        const vector<FactProxy> &pre,
        const std::vector<FactProxy> &pre2) const {
    for (const auto &p : pre){
        if (find_if(begin(pre2), end(pre2),
                    [p] (const FactProxy &p2){
                        return p2 == p;
                    }) == end(pre2)){
            return false;
        }
    }
    return true;
}

void SimulationRelation::compute_label_dominance() {
    op_dominated_by.resize(task->get_num_operators());
    for (auto & elem : op_dominated_by){
        elem.resize(task->get_num_operators(), false);
    }

    vector<FactorID> op_to_leaf(task->get_num_operators(), FactorID::CENTER);
    vector<vector<FactProxy>> op_to_center_pre(task->get_num_operators());

    for (auto op : task_proxy.get_operators()) {
        if (!factoring->is_global_operator(op.get_id())) {
            for (FactorID leaf(0); leaf < factoring->get_num_leaves(); ++leaf) {
                if (factoring->has_eff_on_leaf(op.get_id(), leaf)){
                    op_to_leaf[op.get_id()] = leaf;
                }
            }
            for (auto pre : op.get_preconditions()){
                if (factoring->get_factor(pre.get_variable().get_id()) == FactorID::CENTER){
                    op_to_center_pre[op.get_id()].push_back(pre);
                }
            }
        }
    }

    for (auto op : task_proxy.get_operators()) {
        FactorID factor = op_to_leaf[op.get_id()];
        if (factor == FactorID::CENTER || !factoring->is_fork_leaf(factor)) {
            // skip global and non-fork leaf actions
            continue;
        }

        int cost = op.get_cost();
        const auto &pre = op_to_center_pre[op.get_id()];

        for (int op2_id = op.get_id(); op2_id < task->get_num_operators(); ++op2_id) {
            OperatorProxy op2 = task_proxy.get_operators()[op2_id];
            if(factor != op_to_leaf[op2.get_id()]) {
                // Skip actions of other factors
                continue;
            }

            int cost2 = op2.get_cost();
            const auto &pre2 = op_to_center_pre[op2.get_id()];

            if(cost <= cost2 && center_precondition_dominance(pre, pre2)) {
                op_dominated_by[op2.get_id()][op.get_id()] = true;
            }

            if (cost2 <= cost && center_precondition_dominance(pre2, pre)) {
                op_dominated_by[op.get_id()][op2.get_id()] = true;
            }
        }
    }
}

void SimulationRelation::compute_simulation(FactorID factor, const utils::CountdownTimer &timer) {
    // Init goal respecting
    size_t num_states = leaf_state_space.get_num_states(factor);
    vector<int> goal_distances(num_states, numeric_limits<int>::max());
    priority_queues::BucketQueue<LeafStateHash> open;

    for (LeafStateHash s : leaf_state_space.leaf_goal_states[factor]) {
        goal_distances[s] = 0;
        open.push(0, s);
    }

    while (!open.empty()) {
        if (timer.is_expired()){
            relation[factor].clear();
            return;
        }
        pair<int, LeafStateHash> entry = open.pop();
        LeafStateHash state = entry.second;
        int value = entry.first;
        if (goal_distances[state] < value) {
            continue;
        }
        for (const auto &transition :  leaf_state_space.leaf_state_predecessors[factor][state]) {
            LeafStateHash t = transition.second;

            if (value + 1 < goal_distances[t]) {
                goal_distances[t] = value + 1;
                open.push(value + 1, t);
            }
        }
    }

    relation[factor].resize(num_states);
    for(LeafStateHash i(0); i < num_states; ++i){
        if (timer.is_expired()){
            relation[factor].clear();
            return;
        }
        relation[factor][i].resize(num_states, true);
        if(!leaf_state_space.is_leaf_goal_state(i, factor)){
            for (LeafStateHash j(0); j < num_states; ++j){
                if (leaf_state_space.is_leaf_goal_state(j, factor) || goal_distances[i] > goal_distances[j]){
                    relation[factor][i][j] = false;
                }
            }
        }
    }

    bool changes = true;
    while (changes) {
        changes = false;
        if (timer.is_expired()){
            relation[factor].clear();
            return;
        }
        for (LeafStateHash s(0); s < leaf_state_space.get_num_states(factor); ++s) {
            for (LeafStateHash t(0); t < leaf_state_space.get_num_states(factor); ++t) { //for each pair of states t, s
                if (s != t && simulates(factor, t, s)) {
                    //Check if really t simulates s
                    //for each transition s--l->s':
                    // a) with noop t >= s' and l dominated by noop?
                    // b) exist t--l'-->t', t' >= s' and l dominated by l'?
                    for (const auto &trs : leaf_state_space.leaf_state_successors[factor][s]){
                        LeafStateHash trs_target = trs.second;
                        OperatorID trs_label = trs.first;

                        if(simulates(factor, t, trs_target)) {
                            continue;
                        }
                        bool found = false;
                        for (const auto &trt  : leaf_state_space.leaf_state_successors[factor][t]) {
                            LeafStateHash trt_target = trt.second;
                            OperatorID trt_label = trt.first;

                            if(op_dominated_by[trs_label.hash()][trt_label.hash()] &&
                                    simulates(factor, trt_target, trs_target)) {
                                found = true;
                                break;
                            }
                        }
                        if(!found) {
                            changes = true;
                            remove(factor, t, s);
                        }
                    }
                }
            }
        }
    }
}

size_t SimulationRelation::num_equivalences(FactorID factor) const {
    size_t num = 0;
    vector<bool> counted(relation[factor].size(), false);
    for(LeafStateHash i(0); i < counted.size(); ++i){
        if(!counted[i]){
            for(LeafStateHash j(i + 1); j < relation[factor].size(); ++j){
                if(similar(factor, i, j)){
                    counted [j] = true;
                }
            }
        } else {
            ++num;
        }
    }
    return num;
}

size_t SimulationRelation::num_simulations(FactorID factor, bool ignore_equivalences) const {
    size_t res = 0;
    if (ignore_equivalences){
        vector<bool> counted (relation[factor].size(), false);
        for(LeafStateHash i(0); i < relation[factor].size(); ++i){
            if(!counted[i]){
                for(LeafStateHash j(i+1); j < relation[factor].size(); ++j){
                    if(similar(factor, i, j)){
                        counted[j] = true;
                    }
                }
            }
        }
        for(LeafStateHash i(0); i < relation[factor].size(); ++i){
            if(!counted[i]){
                for(LeafStateHash j(i+1); j < relation[factor].size(); ++j){
                    if(!counted[j]){
                        if(!similar(factor, i, j) && (simulates(factor, i, j) || simulates(factor, j, i))){
                            ++res;
                        }
                    }
                }
            }
        }
    } else {
        for(LeafStateHash i(0); i < relation[factor].size(); ++i){
            for(LeafStateHash j(0); j < relation[factor].size(); ++j){
                if(simulates(factor, i, j)){
                    ++res;
                }
            }
        }
    }
    return res;
}

void SimulationRelation::statistics() const {
    log << "Simulation Relation Finished" << endl;
    if (relation.empty()){
        log << "No relation computed due to timeout." << endl;
    }
    for (FactorID factor(0); factor < relation.size(); ++factor){
        log << "Factor " <<  factor << " has " <<
                num_equivalences(factor) << " equivalences and " <<
                num_simulations(factor, true) << " simulations " << endl;
    }
}

void SimulationRelation::dump(FactorID factor) const {
    if (relation[factor].empty()){
        // is a non-fork leaf or the computation timed out
        return;
    }

    log << "SIMREL:" << endl;

    for(LeafStateHash j(0); j < relation[factor].size(); ++j){
        for(LeafStateHash i(0); i < relation[factor][i].size(); ++i){
            if(simulates(factor, j, i) && i != j){
                if(simulates(factor, i, j)){
                    if (j < i){
                        log << get_name(i, factor) << " <=> " << get_name(j, factor) << endl;
                    }
                } else {
                    log << get_name(i, factor) << " <= " << get_name(j, factor) << endl;
                }
            }
        }
    }

    log << "Reasons: " << endl;
    for (LeafStateHash s(0); s < leaf_state_space.get_num_states(factor); ++s) {
        for (LeafStateHash t(0); t < leaf_state_space.get_num_states(factor); ++t) { //for each pair of states t, s
            if (s != t && simulates(factor, t, s)) {
                for (const auto &trs  : leaf_state_space.leaf_state_successors[factor][s]){
                    LeafStateHash trs_target = trs.second;
                    OperatorID trs_label = trs.first;

                    if(simulates(factor, t, trs_target)) {
                        log << get_name(s, factor) << " -> " << get_name(trs_target, factor) << " is simulated by " <<
                                get_name(t, factor) << " noop " << endl;
                        continue;
                    }
                    for (const auto &trt  : leaf_state_space.leaf_state_successors[factor][t]) {
                        LeafStateHash trt_target = trt.second;
                        OperatorID trt_label = trt.first;

                        if(op_dominated_by[trs_label.hash()][trt_label.hash()] &&
                                simulates(factor, trt_target, trs_target)) {
                            log << get_name(s, factor) << " -> " << get_name(trs_target, factor) << " is simulated by " <<
                                    get_name(t, factor) << " -> " << get_name(trt_target, factor) << endl;
                            break;
                        }
                    }
                }
            }
        }
    }
}

string SimulationRelation::get_name(LeafStateHash id, FactorID factor) const {
    return leaf_state_space.get_name(id, factor);
}
}
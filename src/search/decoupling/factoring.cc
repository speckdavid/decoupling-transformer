#include "factoring.h"

#include "leaf_state.h"
#include "plan_reconstruction.h"
#include "../plugins/plugin.h"
#include "../tasks/root_task.h"
#include "../task_utils/task_properties.h"

using namespace std;

namespace decoupling {
Factoring::Factoring(const plugins::Options &opts) :
    prune_fork_leaf_state_spaces(opts.get<bool>("prune_fork_leaf_state_spaces")),
    num_global_operators(0),
    log(utils::get_log_from_options(opts)),
    factoring_timer(utils::CountdownTimer(opts.get<int>("factoring_time_limit"))),
    task(tasks::g_root_task),
    task_proxy(TaskProxy(*task)),
    min_number_leaves(opts.get<int>("min_number_leaves")),
    max_leaf_size(opts.get<int>("max_leaf_size")) {
    task_properties::verify_no_axioms(task_proxy);
    task_properties::verify_no_conditional_effects(task_proxy);

    if (prune_fork_leaf_state_spaces) {
        log << "Setting prune_fork_leaf_state_spaces=true is not (yet) supported." << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }
}

void Factoring::apply_factoring() {
    // initialize center/leaves
    var_to_factor.resize(task->get_num_variables(), FactorID::CENTER);
    var_to_id_in_factor.resize(task->get_num_variables(), -1);

    // normalize the factoring, to be able to compare different factoring methods
    for (auto &leaf : leaves) {
        sort(leaf.begin(), leaf.end());
    }
    sort(leaves.begin(), leaves.end());

    FactorID factor(0);
    for (const auto &leaf : leaves) {
        int i = 0;
        for (int var : leaf) {
            var_to_factor[var] = factor;
            var_to_id_in_factor[var] = i++;
        }
        ++factor;
    }
    for (VariableProxy var : task_proxy.get_variables()) {
        if (var_to_factor[var.get_id()] == FactorID::CENTER) {
            var_to_id_in_factor[var.get_id()] = static_cast<int>(center.size());
            center.push_back(var.get_id());
        }
    }

    // initialize global/leaf-only operators + interaction graph
    interaction_graph = make_unique<InteractionGraph>(get_num_leaves());
    is_global_operator_.resize(task->get_num_operators(), false);
    leaf_operators.resize(get_num_leaves());
    has_op_leaf_pre.resize(get_num_leaves(), vector<bool>(task->get_num_operators(), false));
    has_op_leaf_eff.resize(get_num_leaves(), vector<bool>(task->get_num_operators(), false));
    for (OperatorProxy op : task_proxy.get_operators()) {
        set<FactorID> pre_factors;
        set<FactorID> eff_factors;

        for (auto pre : op.get_preconditions()) {
            pre_factors.insert(var_to_factor[pre.get_variable().get_id()]);
        }

        for (EffectProxy eff : op.get_effects()) {
            eff_factors.insert(var_to_factor[eff.get_fact().get_variable().get_id()]);
        }
        assert(!eff_factors.empty());

        for (FactorID pre_factor : pre_factors) {
            for (FactorID eff_factor : eff_factors) {
                if (pre_factor != eff_factor) {
                    interaction_graph->add_dependency(pre_factor, eff_factor);
                }
            }
        }
        for (FactorID eff1_factor : eff_factors) {
            for (FactorID eff2_factor : eff_factors) {
                if (eff1_factor != eff2_factor) {
                    interaction_graph->add_dependency(eff1_factor, eff2_factor);
                }
            }
        }

        if (eff_factors.find(FactorID::CENTER) != eff_factors.end()) {
            // effect on center variable
            is_global_operator_[op.get_id()] = true;
        } else if (eff_factors.size() > 1) {
            // effect on more than one factor
            is_global_operator_[op.get_id()] = true;
        } else {
            pre_factors.erase(FactorID::CENTER);
            if (pre_factors.size() > 1) {
                // precondition on more than one leaf
                is_global_operator_[op.get_id()] = true;
            } else if (!pre_factors.empty() && *pre_factors.begin() != *eff_factors.begin()) {
                // precondition on leaf A, but effect on leaf B
                is_global_operator_[op.get_id()] = true;
            }
        }
        assert(is_global_operator_[op.get_id()] || eff_factors.size() == 1);
        if (is_global_operator_[op.get_id()]) {
            num_global_operators++;
        }
        for (FactorID leaf : pre_factors) {
            if (leaf != FactorID::CENTER) {
                has_op_leaf_pre[leaf][op.get_id()] = true;
            }
        }
        for (FactorID leaf : eff_factors) {
            if (leaf != FactorID::CENTER) {
                leaf_operators[leaf].emplace_back(op.get_id());
                has_op_leaf_eff[leaf][op.get_id()] = true;
            }
        }
    }

    // initialize leaf goals
    goals_by_leaf.resize(leaves.size());
    for (auto goal : task_proxy.get_goals()) {
        FactorID leaf = var_to_factor[goal.get_variable().get_id()];
        if (leaf != FactorID::CENTER) {
            goals_by_leaf[leaf].push_back(goal.get_pair());
        }
    }
}

void Factoring::print_factoring() const {
    if (log.is_at_least_verbose()) {
        log << "center factor:" << endl;
        for (int var: center) {
            log << "\t" << task_proxy.get_variables()[var].get_fact(0).get_name() << endl;
        }
        log << "leaf factors:" << endl;
        int i = 0;
        for (const auto &leaf: leaves) {
            log << "leaf " << i++ << ":" << endl;
            for (int var: leaf) {
                log << "\t" << task_proxy.get_variables()[var].get_fact(0).get_name() << endl;
            }
        }
    }
    if (log.is_at_least_normal()) {
        if (interaction_graph->is_fork()) {
            log << "is fork factoring" << endl;
        } else if (interaction_graph->is_ifork()) {
            log << "is inverted-fork factoring" << endl;
        } else if (interaction_graph->is_strict_star()) {
            log << "is strict-star factoring" << endl;
        } else {
            log << "is general factoring" << endl;
        }
    }
}

void Factoring::remove_never_applicable_global_ops(FactorID /*leaf*/) {
    // TODO implement this
}

bool Factoring::does_op_uniquely_fix_lstate(OperatorProxy op, FactorID leaf) const {
    vector<bool> is_var_covered(leaves[leaf].size(), false);
    size_t num_covered_vars = 0;
    for (FactProxy pre : op.get_preconditions()) {
        int var = pre.get_variable().get_id();
        if (get_factor(var) == leaf && !is_var_covered[get_id_in_factor(var)]) {
            is_var_covered[get_id_in_factor(var)] = true;
            num_covered_vars++;
            if (num_covered_vars == leaves[leaf].size()) {
                return true;
            }
        }
    }
    for (EffectProxy eff : op.get_effects()) {
        int var = eff.get_fact().get_variable().get_id();
        if (get_factor(var) == leaf && !is_var_covered[get_id_in_factor(var)]) {
            is_var_covered[get_id_in_factor(var)] = true;
            num_covered_vars++;
            if (num_covered_vars == leaves[leaf].size()) {
                return true;
            }
        }
    }
    return false;
}

bool Factoring::does_op_restrict_leaf(OperatorProxy op, FactorID leaf) const {
    assert(is_global_operator(op.get_id()));
    assert(leaf != FactorID::CENTER && leaf < leaves.size());
    vector<bool> is_op_eff_var(task->get_num_variables(), false);
    for (auto eff : op.get_effects()) {
        is_op_eff_var[eff.get_fact().get_variable().get_id()] = true;
    }
    for (auto op_id : get_leaf_operators(leaf)) {
        if (!is_global_operator(op_id.get_index())) {
            for (auto pre : task_proxy.get_operators()[op_id].get_preconditions()) {
                if (is_op_eff_var[pre.get_variable().get_id()]) {
                    return true;
                }
            }
        }
    }
    return false;
}

bool Factoring::does_op_restrict_leaf(int op_id, int leaf) const {
    return does_op_restrict_leaf(task_proxy.get_operators()[op_id], FactorID(leaf));
}

void Factoring::do_conclusive_leaf_check() {
    is_leaf_conclusive_.resize(leaves.size(), true);
    size_t num_optimizable_leaves = leaves.size();
    for (auto op : task_proxy.get_operators()) {
        if (is_global_operator(op.get_id())) {
            for (FactorID leaf(0); leaf < leaves.size(); ++leaf) {
                if (!is_leaf_conclusive_[leaf]) {
                    continue;
                }
                if (has_pre_or_eff_on_leaf(op.get_id(), leaf)) {
                    // if does_op_uniquely_fix_lstate holds, then after applying op,
                    // there is a unique leaf state reached, which is what we need
                    if (!does_op_uniquely_fix_lstate(op, leaf)) {
                        is_leaf_conclusive_[leaf] = false;
                        num_optimizable_leaves--;
                    }
                } else {
                    if (is_fork_leaf(leaf) && !is_ifork_leaf(leaf)) {
                        // proper fork leafs, i.e. fork leaves with connection to the center, cannot be optimized
                        // this is subsumed by the next check, but cheaper to compute
                        is_leaf_conclusive_[leaf] = false;
                        num_optimizable_leaves--;
                    } else if (does_op_restrict_leaf(op, leaf)) {
                        // the operator restricts the set of reachable leaf states by
                        // en/disabling center preconditions of leaf-only operators
                        is_leaf_conclusive_[leaf] = false;
                        num_optimizable_leaves--;
                    }
                }
            }
            if (num_optimizable_leaves == 0) {
                break;
            }
        }
    }
}

const vector<OperatorID> &Factoring::get_leaf_operators(FactorID leaf) const {
    return leaf_operators[leaf];
}

bool Factoring::has_leaf_goal(FactorID leaf) const {
    assert(leaf != FactorID::CENTER && leaf < leaves.size());
    return !goals_by_leaf[leaf].empty();
}

bool Factoring::is_center_applicable(const State &state, OperatorProxy op) const {
    for (FactProxy pre : op.get_preconditions()) {
        int var = pre.get_variable().get_id();
        int new_id = get_id_in_factor(var);
        if (get_factor(var) == FactorID::CENTER &&
            state[new_id].get_value() != pre.get_value()) {
            return false;
        }
    }
    return true;
}

int Factoring::get_num_effects_on_leaf(OperatorProxy op, FactorID leaf) const {
    int num_effs = 0;
    for (auto eff : op.get_effects()) {
        int var = eff.get_fact().get_variable().get_id();
        if (get_factor(var) == leaf) {
            num_effs++;
        }
    }
    return num_effs;
}

const vector<FactPair> &Factoring::get_leaf_goals(FactorID leaf) const {
    assert(leaf != FactorID::CENTER && leaf < leaves.size());
    return goals_by_leaf[leaf];
}

bool Factoring::is_factoring_possible() const {
    for (auto op : task_proxy.get_operators()) {
        if (static_cast<int>(op.get_effects().size()) < task->get_num_variables()) {
            // there exists a variable that is not affected by all actions
            // => we can construct a mobile factoring with at least one leaf
            return true;
        }
    }
    return false;
}

inline bool is_intersection_empty(const vector<int> &x, const vector<int> &y, int num_vars) {
    if (static_cast<int>(x.size() + y.size()) > num_vars) {
        // this implies a non-empty overlap and is cheap to check
        return false;
    }
    size_t i = 0;
    for (int a : x) {
        for (; i < y.size(); ++i) {
            if (a < y[i]) {
                break;
            } else if (a == y[i]) {
                return false;
            }
        }
    }
    return true;
}

bool Factoring::is_two_leaf_factoring_possible() const {
    int num_vars = task->get_num_variables();
    for (size_t i = 0; i < action_schemas.size(); ++i) {
        const auto &a1 = action_schemas[i];
        for (size_t j = i + 1; j < action_schemas.size(); ++j) {
            const auto &a2 = action_schemas[j];

            if (!is_intersection_empty(a1.eff_vars, a2.eff_vars, num_vars)) {
                // potential leaves overlap
                continue;
            }
            if (!is_intersection_empty(a1.eff_vars, a2.pre_vars, num_vars)) {
                // actions of one leaf are preconditioned by variables in the other leaf
                continue;
            }
            if (!is_intersection_empty(a1.pre_vars, a2.eff_vars, num_vars)) {
                // actions of one leaf are preconditioned by variables in the other leaf
                continue;
            }
            // there exists a mobile factoring with at leaf two leaves: a1.eff_vars and a2.eff_vars
            return true;
        }
    }
    log << "No mobile factoring with at least 2 leaves possible." << endl;
    return false;
}

void Factoring::check_factoring() const {
    // TODO add more sanity checks?
    // check if there are any leaves
    if (leaves.empty()) {
        log << "No factoring with at least " << min_number_leaves << " leaves found." << endl;
        utils::exit_with(utils::ExitCode::SEARCH_UNSOLVED_INCOMPLETE);
    }
    // check that leaves are pairwise disjoint
    vector<bool> is_leaf_var(task->get_num_variables(), false);
    for (const auto &leaf : leaves) {
        for (int var : leaf) {
            if (is_leaf_var[var]) {
                log << "ERROR: invalid factoring, variable "
                    << task_proxy.get_variables()[var].get_fact(0).get_name()
                    << " appears in multiple leaves." << endl;
                utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
            }
            is_leaf_var[var] = true;
        }
    }
}

void Factoring::compute_factoring() {
    if (!is_factoring_possible()) {
        utils::exit_with(utils::ExitCode::SEARCH_UNSOLVED_INCOMPLETE);
    } else if (min_number_leaves > 1) {
        compute_action_schemas();
        if (!is_two_leaf_factoring_possible()) {
            utils::exit_with(utils::ExitCode::SEARCH_UNSOLVED_INCOMPLETE);
        }
    }
    compute_factoring_();
    log << "Factoring time: " << factoring_timer.get_elapsed_time() << endl;
    check_factoring();
    log << "Number leaf factors: " << leaves.size() << endl;
    apply_factoring();

    print_factoring();
    leaf_state_space = make_unique<LeafStateSpace>(shared_from_this(),
                                                   task,
                                                   log,
                                                   false,
                                                   prune_fork_leaf_state_spaces);
    save_memory();
}

void Factoring::save_memory() {
    // TODO save more memory, what about the leaf state space?
    vector<ActionSchema>().swap(action_schemas);
    vector<vector<int>>().swap(var_to_affecting_op);
}

bool Factoring::check_timeout() const {
    if (factoring_timer.is_expired()) {
        return false;
    }
    return true;
}

void Factoring::compute_action_schemas() {
    if (action_schemas.empty()) {
        OperatorsProxy operators = task_proxy.get_operators();
        assert(!operators.empty());
        utils::HashMap<vector<int>, utils::HashMap<vector<int>, size_t>> scheme_loockup;
        for (OperatorProxy op : operators) {
            vector<int> pre_vars;
            for (FactProxy pre : op.get_preconditions()) {
                pre_vars.push_back(pre.get_variable().get_id());
            }
            sort(pre_vars.begin(), pre_vars.end());

            vector<int> eff_vars;
            for (EffectProxy eff : op.get_effects()) {
                eff_vars.push_back(eff.get_fact().get_variable().get_id());
            }
            sort(eff_vars.begin(), eff_vars.end());

            auto it = scheme_loockup.find(pre_vars);
            if (it == scheme_loockup.end() || it->second.find(eff_vars) == it->second.end()) {
                scheme_loockup[pre_vars][eff_vars] = action_schemas.size();
                action_schemas.emplace_back(1, pre_vars, eff_vars);
            } else {
                action_schemas[it->second[eff_vars]].inc_num_actions();
            }
        }
    }
}

void Factoring::compute_var_to_ops_map() {
    if (var_to_affecting_op.empty()) {
        var_to_affecting_op.resize(task->get_num_variables());
        for (OperatorProxy op : task_proxy.get_operators()) {
            for (EffectProxy eff : op.get_effects()) {
                var_to_affecting_op[eff.get_fact().get_variable().get_id()].push_back(op.get_id());
            }
        }
    }
}

int Factoring::get_factor(int var) const {
    return var_to_factor[var];
}

int Factoring::get_id_in_factor(int var) const {
    return var_to_id_in_factor[var];
}

bool Factoring::is_center_variable(int var) const {
    return var_to_factor[var] == FactorID::CENTER;
}

bool Factoring::is_leaf_variable(int var) const {
    return var_to_factor[var] != FactorID::CENTER;
}

bool Factoring::is_leaf_only_operator(int operator_id) const {
    return !is_global_operator_[operator_id];
}

bool Factoring::is_global_operator(int operator_id) const {
    return is_global_operator_[operator_id];
}

int Factoring::get_num_leaves() const {
    return leaves.size();
}

int Factoring::get_num_leaf_variables(int leaf) const {
    assert(leaf < (int)leaves.size());
    return leaves.at(leaf).size();
}

int Factoring::get_num_leaf_states(int leaf_) const {
    FactorID leaf(leaf_);
    assert(leaf < leaves.size());
    return leaf_state_space->get_num_states(leaf);
}

int Factoring::get_num_all_leaf_states() const {
    int res = 0;
    for (FactorID leaf(0); leaf < get_num_leaves(); ++leaf) {
        res += get_num_leaf_states(leaf);
    }
    return res;
}

int Factoring::get_num_all_goal_leaf_states() const {
    int res = 0;
    for (FactorID leaf(0); leaf < get_num_leaves(); ++leaf) {
        res += get_goal_leaf_states(leaf).size();
    }
    return res;
}

int Factoring::get_num_global_operators() const {
    return num_global_operators;
}

bool Factoring::is_fork_leaf(FactorID leaf) const {
    return interaction_graph->is_fork_leaf(leaf);
}

bool Factoring::is_ifork_leaf(FactorID leaf) const {
    return interaction_graph->is_ifork_leaf(leaf);
}

bool Factoring::is_fork_leaf(int leaf_) const {
    assert(leaf_ >= 0 && leaf_ < get_num_leaves() && leaf_ != FactorID::CENTER);
    FactorID leaf(leaf_);
    return interaction_graph->is_fork_leaf(leaf);
}

bool Factoring::is_ifork_leaf(int leaf_) const {
    assert(leaf_ >= 0 && leaf_ < get_num_leaves() && leaf_ != FactorID::CENTER);
    FactorID leaf(leaf_);
    return interaction_graph->is_ifork_leaf(leaf);
}

bool Factoring::is_fork_factoring() const {
    for (FactorID leaf(0); leaf < get_num_leaves(); ++leaf) {
        if (!is_fork_leaf(leaf)) {
            return false;
        }
    }
    return true;
}

bool Factoring::has_pre_on_leaf(OperatorID op_id, FactorID leaf) const {
    assert(leaf != FactorID::CENTER && leaf < leaves.size());
    return has_op_leaf_pre[leaf][op_id.get_index()];
}

bool Factoring::has_pre_on_leaf(int op_id, int leaf) const {
    assert(leaf >= 0 && leaf != FactorID::CENTER && leaf < static_cast<int>(leaves.size()));
    return has_op_leaf_pre[leaf][op_id];
}

bool Factoring::has_eff_on_leaf(OperatorID op_id, FactorID leaf) const {
    assert(leaf != FactorID::CENTER && leaf < leaves.size());
    return has_op_leaf_eff[leaf][op_id.get_index()];
}

bool Factoring::has_eff_on_leaf(int op_id, int leaf) const {
    assert(leaf >= 0 && leaf != FactorID::CENTER && leaf < static_cast<int>(leaves.size()));
    return has_op_leaf_eff[leaf][op_id];
}

bool Factoring::has_pre_or_eff_on_leaf(OperatorID op_id, FactorID leaf) const {
    return has_pre_on_leaf(op_id, leaf) || has_eff_on_leaf(op_id, leaf);
}

bool Factoring::has_pre_or_eff_on_leaf(int op_id, int leaf) const {
    return has_pre_on_leaf(op_id, leaf) || has_eff_on_leaf(op_id, leaf);
}

const vector<int> &Factoring::get_center() const {
    return center;
}

const vector<vector<int>> &Factoring::get_leaves() const {
    return leaves;
}

const vector<int> &Factoring::get_leaf(int leaf) const {
    assert(leaf >= 0 && leaf != FactorID::CENTER && leaf < static_cast<int>(leaves.size()));
    return leaves[leaf];
}

int Factoring::get_initial_leaf_state(int /*leaf*/) const {
    return 0; // by convention in the LeafStateRegistry
}

vector<FactPair> Factoring::get_leaf_state_values(int leaf_, int leaf_state) const {
    assert(leaf_ >= 0 && leaf_ != FactorID::CENTER && leaf_ < static_cast<int>(leaves.size()));
    assert(leaf_state >= 0 && (size_t)leaf_state < LeafStateHash::MAX);

    FactorID leaf(leaf_);
    LeafState lstate(leaf_state_space->get_leaf_state(LeafStateHash(leaf_state), leaf));

    vector<FactPair> res;
    for (int var = 0; var < task->get_num_variables(); ++var) {
        if (leaf_ == get_factor(var)) {
            res.emplace_back(var, lstate[var]);
        }
    }
    assert((int)res.size() == get_num_leaf_variables(leaf));
    return res;
}

bool Factoring::is_reachable_condition(const vector<FactPair> &partial_state) {
    unordered_map<int, vector<FactPair>> leaf_to_partial_state;
    for (const FactPair &fact : partial_state) {
        if (is_center_variable(fact.var)) {
            continue;
        }
        int leaf = get_factor(fact.var);
        leaf_to_partial_state[leaf].push_back(fact);
    }

    for (const auto & [l, facts] : leaf_to_partial_state) {
        FactorID leaf(l);
        bool found_state = false;
        for (LeafStateHash id(0); id < leaf_state_space->get_num_states(leaf); ++id) {
            LeafState lstate(leaf_state_space->get_leaf_state(id, leaf));
            bool is_model = true;
            for (const FactPair &fact : facts) {
                if (lstate[fact.var] != fact.value) {
                    is_model = false;
                    break;
                }
            }
            if (is_model) {
                found_state = true;
                break;
            }
        }
        if (!found_state) {
            return false;
        }
    }

    return true;
}

vector<int> Factoring::get_valid_leaf_states(int leaf_, const vector<FactPair> &partial_state) {
    assert(leaf_ >= 0 && leaf_ != FactorID::CENTER && leaf_ < static_cast<int>(leaves.size()));

    FactorID leaf(leaf_);

    vector<int> satisfying_lstates;
    for (LeafStateHash id(0); id < leaf_state_space->get_num_states(leaf); ++id) {
        LeafState lstate(leaf_state_space->get_leaf_state(id, leaf));
        bool is_model = true;
        for (const FactPair &fact : partial_state) {
            if (leaf_ == get_factor(fact.var) && lstate[fact.var] != fact.value) {
                is_model = false;
                break;
            }
        }
        if (is_model) {
            satisfying_lstates.push_back(static_cast<int>(id));
        }
    }
    return satisfying_lstates;
}

const vector<LeafStateHash> &Factoring::get_goal_leaf_states(int leaf_) const {
    assert(leaf_ >= 0 && leaf_ != FactorID::CENTER && leaf_ < static_cast<int>(leaves.size()));
    FactorID leaf(leaf_);
    return leaf_state_space->get_leaf_goal_states(leaf);
}

bool Factoring::is_conclusive_leaf(FactorID leaf) {
    assert(leaf != FactorID::CENTER && leaf < leaves.size());
    if (is_leaf_conclusive_.empty()) {
        do_conclusive_leaf_check();
    }
    assert(leaf < is_leaf_conclusive_.size());
    return is_leaf_conclusive_[leaf];
}

bool Factoring::is_conclusive_leaf(int leaf) {
    assert(leaf >= 0);
    return is_conclusive_leaf(FactorID(leaf));
}

vector<int> Factoring::get_predecessors(int leaf_, int leaf_state, int operator_id) const {
    assert(leaf_ >= 0 && leaf_ != FactorID::CENTER && leaf_ < static_cast<int>(leaves.size()));
    assert(leaf_state >= 0 && leaf_state < get_num_leaf_states(leaf_));

    FactorID leaf(leaf_);

    if (has_eff_on_leaf(OperatorID(operator_id), leaf)) {
        return leaf_state_space->get_predecessors(leaf, LeafStateHash(leaf_state), OperatorID(operator_id));
    }

    if (has_pre_on_leaf(OperatorID(operator_id), leaf)) {
        // precondition, but no effect, check if op is applicable in leaf_state, if not => no predecessor
        if (!leaf_state_space->is_applicable(LeafStateHash(leaf_state), leaf, task_proxy.get_operators()[operator_id])) {
            return vector<int>{};
        }
    }

    // no precondition or effect => the only predecessor is the state itself
    return vector<int>{leaf_state};
}

void Factoring::add_leaf_facts_to_state(vector<int> &state, int leaf_, int leaf_state_) const {
    FactorID leaf(leaf_);
    LeafStateHash id(leaf_state_);
    assert(leaf != FactorID::CENTER && leaf < leaves.size());
    assert(id != LeafStateHash::MAX && id < leaf_state_space->get_num_states(leaf));

    LeafState lstate = leaf_state_space->get_leaf_state(id, leaf);
    for (int var : leaves[leaf]) {
        state[var] = lstate[var];
    }
}

string Factoring::get_leaf_name(int leaf_) const {
    FactorID leaf(leaf_);
    assert(leaf < leaves.size());
    return to_string(leaf);
}

string Factoring::get_leaf_state_name(int leaf_, int leaf_state) const {
    FactorID leaf(leaf_);
    assert(leaf < leaves.size());
    assert(leaf_state < get_num_leaf_states(leaf));
    return leaf_state_space->get_name(LeafStateHash(leaf_state), leaf);
}

void Factoring::insert_leaf_paths(vector<OperatorID> &path,
                                  vector<State> &states,
                                  const shared_ptr<AbstractTask> &original_root_task) const {
    PathPrices::insert_leaf_actions(*original_root_task, *this, *leaf_state_space, path, states);
}

void Factoring::add_options_to_feature(plugins::Feature &feature) {
    utils::add_log_options_to_feature(feature);
    feature.add_option<int>("min_number_leaves",
                            "The minimum number of leaf factors.",
                            "2"
                            );
    feature.add_option<int>("max_leaf_size",
                            "Maximum domain-size product of variables in a leaf.",
                            "1000000"
                            );
    feature.add_option<int>("factoring_time_limit",
                            "Time limit for computing the factoring.",
                            "30"
                            );
    feature.add_option<bool>("prune_fork_leaf_state_spaces",
                             "Run simulation-based pruning in fork leaves to reduce their state space, not supported, yet.",
                             "false"
                             );
}

static class EvaluatorCategoryPlugin : public plugins::TypedCategoryPlugin<Factoring> {
public:
    EvaluatorCategoryPlugin() : TypedCategoryPlugin("Factoring") {
        document_synopsis(
            "A factoring that can be used for decoupled search.");
    }
}
_category_plugin;
}

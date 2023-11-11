#include "factoring.h"

#include "leaf_state.h"

#include "../plugins/plugin.h"

#include "../tasks/root_task.h"

using namespace std;

namespace decoupling {
vector<Factoring::ActionSchema> Factoring::action_schemas;

vector<set<int>> Factoring::var_to_affecting_op;

Factoring::Factoring(const plugins::Options &opts) :
    ignore_invertible_root_leaves(opts.get<bool>("ignore_invertible_root_leaves")),
    prune_fork_leaf_state_spaces(opts.get<bool>("prune_fork_leaf_state_spaces")),
    log(utils::get_log_from_options(opts)),
    factoring_timer(utils::CountdownTimer(opts.get<int>("factoring_time_limit"))),
    task(tasks::g_root_task),
    task_proxy(TaskProxy(*task)),
    num_global_operators(0),
    min_number_leaves(opts.get<int>("min_number_leaves")),
    max_leaf_size(opts.get<int>("max_leaf_size")) {
}

void Factoring::apply_factoring() {
    // initialize center/leaves
    var_to_factor.resize(task->get_num_variables(), FactorID::CENTER);
    var_to_id_in_factor.resize(task->get_num_variables(), -1);
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
    is_global_operator_.resize(task->get_num_operators(), false);
    interaction_graph = make_unique<InteractionGraph>(leaves.size());
    leaf_operators.resize(get_num_leaves());
    for (OperatorProxy op : task_proxy.get_operators()){

        set<FactorID> pre_factors;
        set<FactorID> eff_factors;

        for (auto pre : op.get_preconditions()){
            pre_factors.insert(var_to_factor[pre.get_variable().get_id()]);
        }

        for (EffectProxy eff : op.get_effects()){
            eff_factors.insert(var_to_factor[eff.get_fact().get_variable().get_id()]);
        }
        assert(!eff_factors.empty());

        for (FactorID pre_factor : pre_factors){
            for (FactorID eff_factor : eff_factors){
                if (pre_factor != eff_factor){
                    interaction_graph->add_dependency(pre_factor, eff_factor);
                }
            }
        }
        for (FactorID eff1_factor : eff_factors){
            for (FactorID eff2_factor : eff_factors){
                if (eff1_factor != eff2_factor){
                    interaction_graph->add_dependency(eff1_factor, eff2_factor);
                }
            }
        }

        if (eff_factors.find(FactorID::CENTER) != eff_factors.end()){
            // effect on center variable
            is_global_operator_[op.get_id()] = true;
        } else if (eff_factors.size() > 1){
            // effect on more than one factor
            is_global_operator_[op.get_id()] = true;
        } else {
            pre_factors.erase(FactorID::CENTER);
            if (pre_factors.size() > 1){
                // precondition on more than one leaf
                is_global_operator_[op.get_id()] = true;
            } else if (!pre_factors.empty() && *pre_factors.begin() != *eff_factors.begin()){
                // precondition on leaf A, but effect on leaf B
                is_global_operator_[op.get_id()] = true;
            }
        }
        assert(is_global_operator_[op.get_id()] || eff_factors.size() == 1);
        if (is_global_operator_[op.get_id()]) {
            num_global_operators++;
        }
        for (FactorID leaf : eff_factors){
            if (leaf != FactorID::CENTER) {
                leaf_operators[leaf].emplace_back(op.get_id());
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
    if (log.is_at_least_normal()) {
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
        if (interaction_graph->is_fork()){
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

const vector<OperatorID> &Factoring::get_leaf_operators(FactorID leaf) const {
    return leaf_operators[leaf];
}

bool Factoring::has_leaf_goal(FactorID leaf) const {
    assert(leaf != FactorID::CENTER);
    return !goals_by_leaf[leaf].empty();
}

const std::vector<FactPair> &Factoring::get_leaf_goals(FactorID leaf) const {
    assert(leaf != FactorID::CENTER);
    return goals_by_leaf[leaf];
}

bool Factoring::is_factoring_possible() const {
    // if there exists a variable that is affected by all actions, then
    // no mobile factoring can exist.
    vector<int> op_count(task->get_num_variables(), 0);
    vector<bool> var_not_affected_by_some_op(task->get_num_variables(), false);
    int num_vars_not_affected_by_some_op = 0;
    for (int i = 0; i < task->get_num_operators(); ++i) {
        for (const EffectProxy &eff : task_proxy.get_operators()[i].get_effects()) {
            int eff_var = eff.get_fact().get_variable().get_id();
            if (op_count[eff_var] == i) {
                ++op_count[eff_var];
            } else if (!var_not_affected_by_some_op[eff_var]) {
                var_not_affected_by_some_op[eff_var] = true;
                ++num_vars_not_affected_by_some_op;
                if (num_vars_not_affected_by_some_op == task->get_num_variables()) {
                    // no variable is affected by all actions
                    return true;
                }
            }
        }
    }
    for (int op_c : op_count) {
        if (op_c == task->get_num_operators()) {
            log << "No mobile factoring possible." << endl;
            return false;
        }
    }
    return true;
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
        log << "ERROR: no factoring with at least " << min_number_leaves << " leaves found." << endl;
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
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
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    } else if (min_number_leaves > 1) {
        compute_action_schemas();
        if (!is_two_leaf_factoring_possible()) {
            utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
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
                                                   ignore_invertible_root_leaves,
                                                   prune_fork_leaf_state_spaces);
}

bool Factoring::check_timeout() const {
    if (factoring_timer.is_expired()) {
        return false;
    }
    return true;
}

void Factoring::compute_action_schemas() {
    if (action_schemas.empty()) {
        OperatorsProxy operators = TaskProxy(*task).get_operators();
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

            if (scheme_loockup.find(pre_vars) == scheme_loockup.end() ||
                scheme_loockup[pre_vars].find(eff_vars) == scheme_loockup[pre_vars].end()) {
                scheme_loockup[pre_vars][eff_vars] = action_schemas.size();
                action_schemas.emplace_back(1, pre_vars, eff_vars);
            } else {
                action_schemas[scheme_loockup[pre_vars][eff_vars]].incr_num_action();
            }
        }
    }
}

void Factoring::compute_var_to_ops_map() {
    if (var_to_affecting_op.empty()) {
        var_to_affecting_op = vector<set<int>>(task_proxy.get_variables().size(), set<int>());
        for (OperatorProxy op : task_proxy.get_operators()) {
            for (EffectProxy eff : op.get_effects()) {
                var_to_affecting_op[eff.get_fact().get_variable().get_id()].insert(op.get_id());
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
    if (is_global_operator_[operator_id]){
        return false;
    }
    set<FactorID> affected_leaves;
    for (EffectProxy eff : task_proxy.get_operators()[operator_id].get_effects()) {
        FactorID var_factor = var_to_factor[eff.get_fact().get_variable().get_id()];
        if (var_factor == FactorID::CENTER) {
            return false;
        } else {
            affected_leaves.insert(var_factor);
            if (affected_leaves.size() > 1){
                return false;
            }
        }
    }
    return true;
}

bool Factoring::is_global_operator(int operator_id) const {
    return is_global_operator_[operator_id];
}

int Factoring::get_num_leaves() const {
    return leaves.size();
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

const vector<int> &Factoring::get_center() const {
    return center;
}

const vector<vector<int>> &Factoring::get_leaves() const {
    return leaves;
}

const vector<int> &Factoring::get_leaf(int leaf_) const {
    FactorID leaf(leaf_);
    assert(leaf < leaves.size());
    return leaves.at(leaf);
}

int Factoring::get_initial_leaf_state(int /*leaf*/) const {
    return 0; // by convention in the LeafStateRegistry
}

vector<int> Factoring::get_valid_precondition_leaf_states(int leaf_, int op_id) const {
    FactorID leaf(leaf_);
    assert(leaf < leaves.size());
    vector<FactProxy> leaf_pre;
    for (FactProxy pre : task_proxy.get_operators()[op_id].get_preconditions()){
        if (var_to_factor[pre.get_variable().get_id()] == leaf){
            leaf_pre.push_back(pre);
        }
    }
    assert(!leaf_pre.empty()); // this should only be called if op_id has a precondition on leaf

    vector<int> satisfying_lstates;
    for (LeafStateHash id(0); id < leaf_state_space->get_num_states(leaf); ++id){
        LeafState lstate(leaf_state_space->get_leaf_state(id, leaf));
        bool applicable = true;
        for (FactProxy pre : leaf_pre){
            if (lstate[pre.get_variable()] != pre.get_value()){
                applicable = false;
                break;
            }
        }
        if (applicable){
            satisfying_lstates.push_back(static_cast<int>(id));
        }
    }

    return satisfying_lstates;
}

const vector<LeafStateHash> &Factoring::get_goal_leaf_states(int leaf_) const {
    FactorID leaf(leaf_);
    assert(leaf < leaves.size());
    return leaf_state_space->get_leaf_goal_states(FactorID(leaf));
}

vector<int> Factoring::get_predecessors(int leaf_, int leaf_state, int operator_id) const {
    FactorID leaf(leaf_);
    assert(leaf < leaves.size());
    assert((int)leaf_state < get_num_leaf_states(leaf));

    // HACK: Check of operators has an effect on the leaf
    // Make this more efficeint and less hacky!
    for (int eff_id = 0; eff_id < task->get_num_operator_effects(operator_id, false); ++eff_id) {
        int eff_var = task->get_operator_effect(operator_id, eff_id, false).var;
        if (get_factor(eff_var) == leaf_) {
            return leaf_state_space->get_predecessors(leaf, LeafStateHash(leaf_state), OperatorID(operator_id));
        }
    }
    
    // IMPORTANT: What about preconditions?
    // If it has no effect on the current leaf, we just return the same state
    return  move(vector<int>{leaf_state});
}

string Factoring::get_leaf_name(int leaf_) const {
    FactorID leaf(leaf_);
    assert(leaf < leaves.size());
    return to_string(leaf);
}

string Factoring::get_leaf_state_name(int leaf_, int leaf_state) const {
    FactorID leaf(leaf_);
    assert(leaf < leaves.size());
    assert((int)leaf_state < get_num_leaf_states(leaf));
    return leaf_state_space->get_name(LeafStateHash(leaf_state), leaf);
}

void Factoring::add_options_to_feature(plugins::Feature &feature) {
    utils::add_log_options_to_feature(feature);
    feature.add_option<int>("min_number_leaves",
                            "maximum number of leaves",
                            "2");
    feature.add_option<int>("max_leaf_size",
                            "maximum domain size product of variables in a leaf",
                            "infinity");
    feature.add_option<int>("factoring_time_limit",
                            "timeout for computing the factoring",
                            "infinity"
                            );
    feature.add_option<bool>("ignore_invertible_root_leaves",
                            "root leaves with invertible state space will be removed",
                            "false"
                            );
    feature.add_option<bool>("prune_fork_leaf_state_spaces",
                             "run simulation-based pruning in fork leaves to reduce their state space",
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

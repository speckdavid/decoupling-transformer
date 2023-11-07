#include "factoring.h"

#include "../plugins/plugin.h"

#include "../tasks/root_task.h"

using namespace std;

namespace decoupling {
vector<Factoring::ActionSchema> Factoring::action_schemas;

vector<set<int>> Factoring::var_to_affecting_op;

Factoring::Factoring(const plugins::Options &opts) :
    log(utils::get_log_from_options(opts)),
    factoring_timer(utils::CountdownTimer(opts.get<int>("factoring_time_limit"))),
    task(tasks::g_root_task),
    task_proxy(TaskProxy(*task)),
    min_number_leaves(opts.get<int>("min_number_leaves")),
    max_leaf_size(opts.get<int>("max_leaf_size")) {
}

void Factoring::apply_factoring() {
    // TODO implement
    var_to_factor.resize(task->get_num_variables(), FactorID::CENTER);
    FactorID factor(0);
    for (const auto &leaf : leaves){
        for (int var : leaf){
            var_to_factor[var] = factor;
        }
        ++factor;
    }
    for (VariableProxy var : task_proxy.get_variables()){
        if (var_to_factor[var.get_id()] == FactorID::CENTER){
            center.push_back(var.get_id());
        }
    }
}

void Factoring::print_factoring() const {
    if (log.is_at_least_normal()) {
        log << "factoring with " << leaves.size() << " leaves" << endl;
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
}

bool Factoring::is_factoring_possible() const {
    // if there exists a variable that is affected by all actions, then
    // no mobile factoring can exist.
    vector<int> op_count(task->get_num_variables(), 0);
    vector<bool> var_not_affected_by_some_op(task->get_num_variables(), false);
    int num_vars_not_affected_by_some_op = 0;
    for (int i = 0; i < task->get_num_operators(); ++i){
        for (const EffectProxy &eff : task_proxy.get_operators()[i].get_effects()){
            int eff_var = eff.get_fact().get_variable().get_id();
            if (op_count[eff_var] == i){
                ++op_count[eff_var];
            } else if (!var_not_affected_by_some_op[eff_var]) {
                var_not_affected_by_some_op[eff_var] = true;
                ++num_vars_not_affected_by_some_op;
                if (num_vars_not_affected_by_some_op == task->get_num_variables()){
                    // no variable is affected by all actions
                    return true;
                }
            }
        }
    }
    for (int op_c : op_count){
        if (op_c == task->get_num_operators()){
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
    for (int a : x){
        for (; i < y.size(); ++i){
            if (a < y[i]){
                break;
            } else if (a == y[i]){
                return false;
            }
        }
    }
    return true;
}

bool Factoring::is_two_leaf_factoring_possible() const {
    int num_vars = task->get_num_variables();
    for (size_t i = 0; i < action_schemas.size(); ++i){
        const auto &a1 = action_schemas[i];
        for (size_t j = i+1; j < action_schemas.size(); ++j){
            const auto &a2 = action_schemas[j];

            if (!is_intersection_empty(a1.eff_vars, a2.eff_vars, num_vars)){
                // potential leaves overlap
                continue;
            }
            if (!is_intersection_empty(a1.eff_vars, a2.pre_vars, num_vars)){
                // actions of one leaf are preconditioned by variables in the other leaf
                continue;
            }
            if (!is_intersection_empty(a1.pre_vars, a2.eff_vars, num_vars)){
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
    if (leaves.empty()){
        log << "ERROR: no factoring with at least " << min_number_leaves << " leaves found." << endl;
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    }
    // check that leaves are pairwise disjoint
    vector<bool> is_leaf_var(task->get_num_variables(), false);
    for (const auto &leaf : leaves){
        for (int var : leaf){
            if (is_leaf_var[var]){
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
    if (!is_factoring_possible()){
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    } else if (min_number_leaves > 1) {
        compute_action_schemas();
        if (!is_two_leaf_factoring_possible()){
            utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
        }
    }
    compute_factoring_();
    log << "Factoring time: " << factoring_timer.get_elapsed_time() << endl;
    check_factoring();
    log << "Number leaf factors: " << leaves.size() << endl;
    apply_factoring();
    print_factoring();
}

bool Factoring::check_timeout() const {
    if (factoring_timer.is_expired()){
        return false;
    }
    return true;
}

void Factoring::compute_action_schemas() {
    if (action_schemas.empty()){
        OperatorsProxy operators = TaskProxy(*task).get_operators();
        assert(!operators.empty());
        utils::HashMap<std::vector<int>, utils::HashMap<std::vector<int>, size_t> > scheme_loockup;
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
                    scheme_loockup[pre_vars].find(eff_vars) == scheme_loockup[pre_vars].end()){
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
        var_to_affecting_op = vector<set<int> >(task_proxy.get_variables().size(), set<int>());
        for (OperatorProxy op : task_proxy.get_operators()) {
            for (EffectProxy eff : op.get_effects()) {
                var_to_affecting_op[eff.get_fact().get_variable().get_id()].insert(op.get_id());
            }
        }
    }
}


bool Factoring::is_center_variable(int var) const {
    return std::find(center.begin(), center.end(), var) != center.end();
}

bool Factoring::is_leaf_variable(int var) const {
    return !is_center_variable(var);
}

int Factoring::get_num_leaves() const {
    return leaves.size();
}

int Factoring::get_num_leaf_states(size_t l) const {
    assert(l < leaves.size());
    // TODO: fix
    int num_states = 1;
    for (int var : leaves.at(l)) {
        num_states *= task->get_variable_domain_size(var);
    }
    return num_states;
}

vector<int> Factoring::get_center() const {
    return center;
}

std::vector<std::vector<int>> Factoring::get_leaves() const {
    return leaves;
}

vector<int> Factoring::get_leaf(size_t l) const {
    assert(l < leaves.size());
    return leaves.at(l);
}

int Factoring::get_initial_leaf_state(size_t l) const {
    assert(l < leaves.size());
    // TODO: return the leaf state id which corresponds to the initial state
    return 0;
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

#include "factoring.h"

#include "../plugins/plugin.h"

#include "../tasks/root_task.h"

using namespace std;

namespace decoupling {
vector<Factoring::ActionSchema> Factoring::action_schemas;

vector<set<int>> Factoring::var_to_affecting_op;

Factoring::Factoring(const plugins::Options &opts) :
    factoring_timer(utils::CountdownTimer(opts.get<int>("factoring_time_limit"))),
    task(tasks::g_root_task),
    task_proxy(TaskProxy(*task)),
    min_number_leaves(opts.get<int>("min_number_leaves")),
    max_leaf_size(opts.get<int>("max_leaf_size")) {
}

void Factoring::apply_factoring() {
    // TODO implement
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

            if (scheme_loockup.find(pre_vars) == scheme_loockup.end()){
                scheme_loockup[pre_vars][eff_vars] = action_schemas.size();
                action_schemas.emplace_back(1, pre_vars, eff_vars);
            } else if (scheme_loockup[pre_vars].find(eff_vars) == scheme_loockup[pre_vars].end()){
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
        TaskProxy task_proxy(*task);
        OperatorsProxy operators = task_proxy.get_operators();
        var_to_affecting_op = vector<set<int> >(task_proxy.get_variables().size(), set<int>());
        for (OperatorProxy op : operators) {
            for (EffectProxy eff : op.get_effects()) {
                var_to_affecting_op[eff.get_fact().get_variable().get_id()].insert(op.get_id());
            }
        }
    }
}

void Factoring::add_options_to_feature(plugins::Feature &feature) {
    feature.add_option<int>("min_number_leaves",
                             "maximum number of leaves",
                             "2");
    feature.add_option<int>("max_leaf_size",
                             "maximum domain size product of variables in a leaf",
                             to_string(numeric_limits<int>::max() - 1));
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

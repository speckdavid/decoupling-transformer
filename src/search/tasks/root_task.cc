#include "root_task.h"

#include "../state_registry.h"

#include "../plugins/plugin.h"
#include "../utils/collections.h"
#include "../utils/timer.h"
#include "../task_utils/task_properties.h"

#include <algorithm>
#include <cassert>
#include <memory>
#include <unordered_set>


using namespace std;
using utils::ExitCode;

namespace tasks {
shared_ptr<AbstractTask> g_root_task = nullptr;

static void check_fact(const FactPair &fact, const vector<ExplicitVariable> &variables) {
    if (!utils::in_bounds(fact.var, variables)) {
        cerr << "Invalid variable id: " << fact.var << endl;
        utils::exit_with(ExitCode::SEARCH_INPUT_ERROR);
    }
    if (fact.value < 0 || fact.value >= variables[fact.var].domain_size) {
        cerr << "Invalid value for variable " << fact.var << ": " << fact.value << endl;
        utils::exit_with(ExitCode::SEARCH_INPUT_ERROR);
    }
}

static void check_facts(const vector<FactPair> &facts, const vector<ExplicitVariable> &variables) {
    for (FactPair fact : facts) {
        check_fact(fact, variables);
    }
}

static void check_facts(const ExplicitOperator &action, const vector<ExplicitVariable> &variables) {
    check_facts(action.preconditions, variables);
    for (const ExplicitEffect &eff : action.effects) {
        check_fact(eff.fact, variables);
        check_facts(eff.conditions, variables);
    }
}

static void check_magic(istream &in, const string &magic) {
    string word;
    in >> word;
    if (word != magic) {
        cerr << "Failed to match magic word '" << magic << "'." << endl
             << "Got '" << word << "'." << endl;
        if (magic == "begin_version") {
            cerr << "Possible cause: you are running the planner "
                 << "on a translator output file from " << endl
                 << "an older version." << endl;
        }
        utils::exit_with(ExitCode::SEARCH_INPUT_ERROR);
    }
}

static vector<FactPair> read_facts(istream &in) {
    int count;
    in >> count;
    vector<FactPair> conditions;
    conditions.reserve(count);
    for (int i = 0; i < count; ++i) {
        FactPair condition = FactPair::no_fact;
        in >> condition.var >> condition.value;
        conditions.push_back(condition);
    }
    return conditions;
}

ExplicitVariable::ExplicitVariable(istream &in) {
    check_magic(in, "begin_variable");
    in >> name;
    in >> axiom_layer;
    in >> domain_size;
    in >> ws;
    fact_names.resize(domain_size);
    for (int i = 0; i < domain_size; ++i)
        getline(in, fact_names[i]);
    check_magic(in, "end_variable");
}

ExplicitVariable::ExplicitVariable(const string &name,
                                   vector<string> &&fact_names,
                                   int axiom_layer,
                                   int axiom_default_value) :
    domain_size(fact_names.size()),
    name(name),
    fact_names(move(fact_names)),
    axiom_layer(axiom_layer),
    axiom_default_value(axiom_default_value) {
}

int ExplicitVariable::get_encoding_size() const {
    return domain_size;
}

ExplicitEffect::ExplicitEffect(
    int var, int value, vector<FactPair> &&conditions)
    : fact(var, value), conditions(move(conditions)) {
}

bool ExplicitEffect::operator==(const ExplicitEffect &other) const {
    return fact == other.fact
           && set<FactPair>(conditions.begin(), conditions.end())
           == set<FactPair>(other.conditions.begin(), other.conditions.end());
}

bool ExplicitEffect::operator<(const ExplicitEffect &other) const {
    if (fact < other.fact)
        return true;
    if (other.fact < fact)
        return false;
    return set<FactPair>(conditions.begin(), conditions.end())
           < set<FactPair>(other.conditions.begin(), other.conditions.end());
}

int ExplicitEffect::get_encoding_size() const {
    return 1 + conditions.size();
}


void ExplicitOperator::read_pre_post(istream &in) {
    vector<FactPair> conditions = read_facts(in);
    int var, value_pre, value_post;
    in >> var >> value_pre >> value_post;
    if (value_pre != -1) {
        preconditions.emplace_back(var, value_pre);
    }
    effects.emplace_back(var, value_post, move(conditions));
}

ExplicitOperator::ExplicitOperator(istream &in, bool is_an_axiom, bool use_metric)
    : is_an_axiom(is_an_axiom) {
    if (!is_an_axiom) {
        check_magic(in, "begin_operator");
        in >> ws;
        getline(in, name);
        preconditions = read_facts(in);
        int count;
        in >> count;
        effects.reserve(count);
        for (int i = 0; i < count; ++i) {
            read_pre_post(in);
        }

        int op_cost;
        in >> op_cost;
        cost = use_metric ? op_cost : 1;
        check_magic(in, "end_operator");
    } else {
        name = "<axiom>";
        cost = 0;
        check_magic(in, "begin_rule");
        read_pre_post(in);
        check_magic(in, "end_rule");
    }
    assert(cost >= 0);
}

ExplicitOperator::ExplicitOperator(int cost,
                                   const string &name,
                                   bool is_an_axiom) :
    cost(cost),
    name(name),
    is_an_axiom(is_an_axiom) {}

ExplicitOperator::ExplicitOperator(int cost,
                                   const string &name,
                                   const vector<FactPair> &preconditions) :
        preconditions(preconditions),
        cost(cost),
        name(name),
        is_an_axiom(false) {}

bool ExplicitOperator::operator==(const ExplicitOperator &other) const {
    return cost == other.cost
           && is_an_axiom == other.is_an_axiom
           && set<FactPair>(preconditions.begin(), preconditions.end())
           == set<FactPair>(other.preconditions.begin(), other.preconditions.end())
           && set<ExplicitEffect>(effects.begin(), effects.end())
           == set<ExplicitEffect>(other.effects.begin(), other.effects.end());
}

bool ExplicitOperator::operator<(const ExplicitOperator &other) const {
    if (cost < other.cost)
        return true;
    if (other.cost < cost)
        return false;
    if (is_an_axiom < other.is_an_axiom)
        return true;
    if (other.is_an_axiom < is_an_axiom)
        return false;

    auto preconditions_set = set<FactPair>(preconditions.begin(), preconditions.end());
    auto other_preconditions_set = set<FactPair>(other.preconditions.begin(), other.preconditions.end());
    if (preconditions_set < other_preconditions_set)
        return true;
    if (other_preconditions_set < preconditions_set)
        return false;

    auto effects_set = set<ExplicitEffect>(effects.begin(), effects.end());
    auto other_effects_set = set<ExplicitEffect>(other.effects.begin(), other.effects.end());
    return effects_set < other_effects_set;
}

int ExplicitOperator::get_encoding_size() const {
    int size = 1 + preconditions.size();
    for (const auto &eff : effects)
        size += eff.get_encoding_size();
    return size;
}

static void read_and_verify_version(istream &in) {
    int version;
    check_magic(in, "begin_version");
    in >> version;
    check_magic(in, "end_version");
    if (version != PRE_FILE_VERSION) {
        cerr << "Expected translator output file version " << PRE_FILE_VERSION
             << ", got " << version << "." << endl
             << "Exiting." << endl;
        utils::exit_with(ExitCode::SEARCH_INPUT_ERROR);
    }
}


static bool read_metric(istream &in) {
    bool use_metric;
    check_magic(in, "begin_metric");
    in >> use_metric;
    check_magic(in, "end_metric");
    return use_metric;
}

static vector<ExplicitVariable> read_variables(istream &in) {
    int count;
    in >> count;
    vector<ExplicitVariable> variables;
    variables.reserve(count);
    for (int i = 0; i < count; ++i) {
        variables.emplace_back(in);
    }
    return variables;
}

static vector<vector<set<FactPair>>> read_mutexes(istream &in, const vector<ExplicitVariable> &variables) {
    vector<vector<set<FactPair>>> inconsistent_facts(variables.size());
    for (size_t i = 0; i < variables.size(); ++i)
        inconsistent_facts[i].resize(variables[i].domain_size);

    int num_mutex_groups;
    in >> num_mutex_groups;

    /*
      NOTE: Mutex groups can overlap, in which case the same mutex
      should not be represented multiple times. The current
      representation takes care of that automatically by using sets.
      If we ever change this representation, this is something to be
      aware of.
    */
    for (int i = 0; i < num_mutex_groups; ++i) {
        check_magic(in, "begin_mutex_group");
        int num_facts;
        in >> num_facts;
        vector<FactPair> invariant_group;
        invariant_group.reserve(num_facts);
        for (int j = 0; j < num_facts; ++j) {
            int var;
            int value;
            in >> var >> value;
            invariant_group.emplace_back(var, value);
        }
        check_magic(in, "end_mutex_group");
        for (const FactPair &fact1 : invariant_group) {
            for (const FactPair &fact2 : invariant_group) {
                if (fact1.var != fact2.var) {
                    /* The "different variable" test makes sure we
                       don't mark a fact as mutex with itself
                       (important for correctness) and don't include
                       redundant mutexes (important to conserve
                       memory). Note that the translator (at least
                       with default settings) removes mutex groups
                       that contain *only* redundant mutexes, but it
                       can of course generate mutex groups which lead
                       to *some* redundant mutexes, where some but not
                       all facts talk about the same variable. */
                    inconsistent_facts[fact1.var][fact1.value].insert(fact2);
                }
            }
        }
    }
    return inconsistent_facts;
}

static vector<FactPair> read_goal(istream &in) {
    check_magic(in, "begin_goal");
    vector<FactPair> goals = read_facts(in);
    check_magic(in, "end_goal");
    if (goals.empty()) {
        cerr << "Task has no goal condition!" << endl;
        utils::exit_with(ExitCode::SEARCH_INPUT_ERROR);
    }
    return goals;
}

static vector<ExplicitOperator> read_actions(
    istream &in, bool is_axiom, bool use_metric,
    const vector<ExplicitVariable> &variables) {
    int count;
    in >> count;
    vector<ExplicitOperator> actions;
    actions.reserve(count);
    for (int i = 0; i < count; ++i) {
        actions.emplace_back(in, is_axiom, use_metric);
        check_facts(actions.back(), variables);
    }
    return actions;
}

RootTask::RootTask(istream &in) {
    read_and_verify_version(in);
    bool use_metric = read_metric(in);
    variables = read_variables(in);
    int num_variables = variables.size();

    mutexes = read_mutexes(in, variables);

    initial_state_values.resize(num_variables);
    check_magic(in, "begin_state");
    for (int i = 0; i < num_variables; ++i) {
        in >> initial_state_values[i];
    }
    check_magic(in, "end_state");

    for (int i = 0; i < num_variables; ++i) {
        variables[i].axiom_default_value = initial_state_values[i];
    }

    goals = read_goal(in);
    check_facts(goals, variables);
    operators = read_actions(in, false, use_metric, variables);
    axioms = read_actions(in, true, use_metric, variables);
    /* TODO: We should be stricter here and verify that we
       have reached the end of "in". */

    /*
      HACK: We use a TaskProxy to access g_axiom_evaluators here which assumes
      that this task is completely constructed.
    */
    AxiomEvaluator &axiom_evaluator = g_axiom_evaluators[TaskProxy(*this)];
    axiom_evaluator.evaluate(initial_state_values);
}

const ExplicitVariable &RootTask::get_variable(int var) const {
    assert(utils::in_bounds(var, variables));
    return variables[var];
}

const ExplicitEffect &RootTask::get_effect(
    int op_id, int effect_id, bool is_axiom) const {
    const ExplicitOperator &op = get_operator_or_axiom(op_id, is_axiom);
    assert(utils::in_bounds(effect_id, op.effects));
    return op.effects[effect_id];
}

const ExplicitOperator &RootTask::get_operator_or_axiom(
    int index, bool is_axiom) const {
    if (is_axiom) {
        assert(utils::in_bounds(index, axioms));
        return axioms[index];
    } else {
        assert(utils::in_bounds(index, operators));
        return operators[index];
    }
}

int RootTask::get_num_variables() const {
    return variables.size();
}

string RootTask::get_variable_name(int var) const {
    return get_variable(var).name;
}

int RootTask::get_variable_domain_size(int var) const {
    return get_variable(var).domain_size;
}

int RootTask::get_variable_axiom_layer(int var) const {
    return get_variable(var).axiom_layer;
}

int RootTask::get_variable_default_axiom_value(int var) const {
    return get_variable(var).axiom_default_value;
}

string RootTask::get_fact_name(const FactPair &fact) const {
    assert(utils::in_bounds(fact.value, get_variable(fact.var).fact_names));
    return get_variable(fact.var).fact_names[fact.value];
}

bool RootTask::are_facts_mutex(const FactPair &fact1, const FactPair &fact2) const {
    if (fact1.var == fact2.var) {
        // Same variable: mutex iff different value.
        return fact1.value != fact2.value;
    }
    assert(utils::in_bounds(fact1.var, mutexes));
    assert(utils::in_bounds(fact1.value, mutexes[fact1.var]));
    return bool(mutexes[fact1.var][fact1.value].count(fact2));
}

int RootTask::get_operator_cost(int index, bool is_axiom) const {
    return get_operator_or_axiom(index, is_axiom).cost;
}

string RootTask::get_operator_name(int index, bool is_axiom) const {
    return get_operator_or_axiom(index, is_axiom).name;
}

int RootTask::get_num_operators() const {
    return operators.size();
}

int RootTask::get_num_operator_preconditions(int index, bool is_axiom) const {
    return get_operator_or_axiom(index, is_axiom).preconditions.size();
}

FactPair RootTask::get_operator_precondition(
    int op_index, int fact_index, bool is_axiom) const {
    const ExplicitOperator &op = get_operator_or_axiom(op_index, is_axiom);
    assert(utils::in_bounds(fact_index, op.preconditions));
    return op.preconditions[fact_index];
}

int RootTask::get_num_operator_effects(int op_index, bool is_axiom) const {
    return get_operator_or_axiom(op_index, is_axiom).effects.size();
}

int RootTask::get_num_operator_effect_conditions(
    int op_index, int eff_index, bool is_axiom) const {
    return get_effect(op_index, eff_index, is_axiom).conditions.size();
}

FactPair RootTask::get_operator_effect_condition(
    int op_index, int eff_index, int cond_index, bool is_axiom) const {
    const ExplicitEffect &effect = get_effect(op_index, eff_index, is_axiom);
    assert(utils::in_bounds(cond_index, effect.conditions));
    return effect.conditions[cond_index];
}

FactPair RootTask::get_operator_effect(
    int op_index, int eff_index, bool is_axiom) const {
    return get_effect(op_index, eff_index, is_axiom).fact;
}

int RootTask::convert_operator_index(
    int index, const AbstractTask *ancestor_task) const {
    if (this != ancestor_task) {
        ABORT("Invalid operator ID conversion");
    }
    return index;
}

int RootTask::get_num_axioms() const {
    return axioms.size();
}

int RootTask::get_num_goals() const {
    return goals.size();
}

FactPair RootTask::get_goal_fact(int index) const {
    assert(utils::in_bounds(index, goals));
    return goals[index];
}

vector<int> RootTask::get_initial_state_values() const {
    return initial_state_values;
}

void RootTask::convert_ancestor_state_values(
    vector<int> &, const AbstractTask *ancestor_task) const {
    if (this != ancestor_task) {
        ABORT("Invalid state conversion");
    }
}

int RootTask::get_encoding_size(bool with_mutexes) const {
    int task_size = 0;
    for (const auto &var : variables)
        task_size += var.get_encoding_size();
    task_size += initial_state_values.size();
    task_size += goals.size();
    for (const auto &op : operators)
        task_size += op.get_encoding_size();
    for (const auto &ax : axioms)
        task_size += ax.get_encoding_size();
    if (with_mutexes) {
        for (size_t var = 0; var < mutexes.size(); ++var) {
            for (size_t val = 0; val < mutexes.at(var).size(); ++val)
                task_size += mutexes.at(var).at(val).size();
        }
    }
    return task_size;
}

void RootTask::normalize_task() {
    sort(goals.begin(), goals.end());
    for (auto &op : operators) {
        sort(op.preconditions.begin(), op.preconditions.end());
        for (auto &eff : op.effects) {
            sort(eff.conditions.begin(), eff.conditions.end());
        }
        sort(op.effects.begin(), op.effects.end());
    }
    for (auto &ax : axioms) {
        sort(ax.preconditions.begin(), ax.preconditions.end());
        for (auto &eff : ax.effects) {
            sort(eff.conditions.begin(), eff.conditions.end());
        }
        sort(ax.effects.begin(), ax.effects.end());
    }
    for (auto &mutex_group : mutexes) {
        sort(mutex_group.begin(), mutex_group.end());
    }
}

void read_root_task(istream &in) {
    assert(!g_root_task);
    g_root_task = make_shared<RootTask>(in);
}

class RootTaskFeature : public plugins::TypedFeature<AbstractTask, AbstractTask> {
public:
    RootTaskFeature() : TypedFeature("no_transform") {
    }

    virtual shared_ptr<AbstractTask> create_component(const plugins::Options &, const utils::Context &) const override {
        return g_root_task;
    }
};

static plugins::FeaturePlugin<RootTaskFeature> _plugin;
}

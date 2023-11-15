#include "undecoupled_task.h"

#include "../plugins/plugin.h"

using namespace std;

namespace tasks {
UndecoupledTask::UndecoupledTask(const std::shared_ptr<AbstractTask> &parent) :
    DelegatingTask(parent),
    decoupled_task(dynamic_pointer_cast<DecoupledRootTask>(parent)) {
    if (!decoupled_task) {
        cerr << "Root task is no decoupled task!" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }
    original_task = decoupled_task->get_original_root_task();
}

int UndecoupledTask::get_num_variables() const {
    return original_task->get_num_variables();
}

std::string UndecoupledTask::get_variable_name(int var) const {
    return original_task->get_variable_name(var);
}

int UndecoupledTask::get_variable_domain_size(int var) const {
    return original_task->get_variable_domain_size(var);
}

int UndecoupledTask::get_variable_axiom_layer(int var) const {
    return original_task->get_variable_axiom_layer(var);
}

int UndecoupledTask::get_variable_default_axiom_value(int var) const {
    return original_task->get_variable_default_axiom_value(var);
}

std::string UndecoupledTask::get_fact_name(const FactPair &fact) const {
    return original_task->get_fact_name(fact);
}

bool UndecoupledTask::are_facts_mutex(
    const FactPair &fact1, const FactPair &fact2) const {
    return original_task->are_facts_mutex(fact1, fact2);
}

int UndecoupledTask::get_operator_cost(int index, bool is_axiom) const {
    return original_task->get_operator_cost(index, is_axiom);
}


std::string UndecoupledTask::get_operator_name(int index, bool is_axiom) const {
    return original_task->get_operator_name(index, is_axiom);
}

int UndecoupledTask::get_num_operators() const {
    return original_task->get_num_operators();
}

int UndecoupledTask::get_num_operator_preconditions(int index, bool is_axiom) const {
    return original_task->get_num_operator_preconditions(index, is_axiom);
}

FactPair UndecoupledTask::get_operator_precondition(
    int op_index, int fact_index, bool is_axiom) const {
    return original_task->get_operator_precondition(op_index, fact_index, is_axiom);
}

int UndecoupledTask::get_num_operator_effects(int op_index, bool is_axiom) const {
    return original_task->get_num_operator_effects(op_index, is_axiom);
}

int UndecoupledTask::get_num_operator_effect_conditions(
    int op_index, int eff_index, bool is_axiom) const {
    return original_task->get_num_operator_effect_conditions(op_index, eff_index, is_axiom);
}

FactPair UndecoupledTask::get_operator_effect_condition(
    int op_index, int eff_index, int cond_index, bool is_axiom) const {
    return original_task->get_operator_effect_condition(op_index, eff_index, cond_index, is_axiom);
}

FactPair UndecoupledTask::get_operator_effect(
    int op_index, int eff_index, bool is_axiom) const {
    return original_task->get_operator_effect(op_index, eff_index, is_axiom);
}

// TODO: fix this
int UndecoupledTask::convert_operator_index_to_parent(int index) const {
    return index;
}

int UndecoupledTask::get_num_axioms() const {
    return original_task->get_num_axioms();
}

int UndecoupledTask::get_num_goals() const {
    return original_task->get_num_goals();
}

FactPair UndecoupledTask::get_goal_fact(int index) const {
    return original_task->get_goal_fact(index);
}

std::vector<int> UndecoupledTask::get_initial_state_values() const {
    return original_task->get_initial_state_values();
}

void UndecoupledTask::convert_state_values_from_parent(std::vector<int> &) const {}

class UndecoupledTaskFeature : public plugins::TypedFeature<AbstractTask, UndecoupledTask> {
public:
    UndecoupledTaskFeature() : TypedFeature("undecoupled") {
        document_title("Undecoupled task");
        document_synopsis(
            "Undecouples a decoupled task to the original one.");
    }

    virtual shared_ptr<UndecoupledTask> create_component(const plugins::Options &/*options*/, const utils::Context &) const override {
        return make_shared<UndecoupledTask>(g_root_task);
    }
};

static plugins::FeaturePlugin<UndecoupledTaskFeature> _plugin;

}

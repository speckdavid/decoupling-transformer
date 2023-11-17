#ifndef TASKS_UNDECOUPLED_TASK_H
#define TASKS_UNDECOUPLED_TASK_H

#include "decoupled_root_task.h"
#include "cost_adapted_task.h"

namespace plugins {
class Options;
}

namespace tasks {

class UndecoupledTask : public CostAdaptedTask {
    std::shared_ptr<DecoupledRootTask> decoupled_task;
    std::shared_ptr<AbstractTask> original_task;

public:
    UndecoupledTask(const std::shared_ptr<AbstractTask> &parent, OperatorCost cost_type);
    virtual ~UndecoupledTask() override = default;

    virtual int get_num_variables() const override;
    virtual std::string get_variable_name(int var) const override;
    virtual int get_variable_domain_size(int var) const override;
    virtual int get_variable_axiom_layer(int var) const override;
    virtual int get_variable_default_axiom_value(int var) const override;
    virtual std::string get_fact_name(const FactPair &fact) const override;
    virtual bool are_facts_mutex(
        const FactPair &fact1, const FactPair &fact2) const override;

    virtual int get_operator_cost(int index, bool is_axiom) const override;
    virtual std::string get_operator_name(int index, bool is_axiom) const override;
    virtual int get_num_operators() const override;
    virtual int get_num_operator_preconditions(int index, bool is_axiom) const override;
    virtual FactPair get_operator_precondition(
        int op_index, int fact_index, bool is_axiom) const override;
    virtual int get_num_operator_effects(int op_index, bool is_axiom) const override;
    virtual int get_num_operator_effect_conditions(
        int op_index, int eff_index, bool is_axiom) const override;
    virtual FactPair get_operator_effect_condition(
        int op_index, int eff_index, int cond_index, bool is_axiom) const override;
    virtual FactPair get_operator_effect(
        int op_index, int eff_index, bool is_axiom) const override;

    virtual int convert_operator_index_to_parent(int index) const override;

    virtual int get_num_axioms() const override;

    virtual int get_num_goals() const override;
    virtual FactPair get_goal_fact(int index) const override;

    virtual std::vector<int> get_initial_state_values() const override;

    virtual void convert_state_values_from_parent(std::vector<int> &) const override;

    void get_sampled_states(const State& dec_state, int samples, std::vector<State>& sampled_states) const;
};
}

#endif
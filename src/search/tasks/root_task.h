#ifndef TASKS_ROOT_TASK_H
#define TASKS_ROOT_TASK_H

#include "../abstract_task.h"

#include <vector>
#include <set>

namespace tasks {
class DecoupledRootTask;

static const int PRE_FILE_VERSION = 3;
extern std::shared_ptr<AbstractTask> g_root_task;
extern void read_root_task(std::istream &in);

struct ExplicitVariable {
    int domain_size;
    std::string name;
    std::vector<std::string> fact_names;
    int axiom_layer;
    int axiom_default_value;

    explicit ExplicitVariable(std::istream &in);
    ExplicitVariable(const std::string &name,
                     std::vector<std::string> &&fact_names,
                     int axiom_layer,
                     int axiom_default_value = 0);

    int get_encoding_size() const;
};


struct ExplicitEffect {
    FactPair fact;
    std::vector<FactPair> conditions;

    ExplicitEffect(int var, int value, std::vector<FactPair> &&conditions);

    bool operator==(const ExplicitEffect &other) const;
    bool operator<(const ExplicitEffect &other) const;

    int get_encoding_size() const;
};


struct ExplicitOperator {
    std::vector<FactPair> preconditions;
    std::vector<ExplicitEffect> effects;
    int cost;
    std::string name;
    bool is_an_axiom;

    void read_pre_post(std::istream &in);
    ExplicitOperator(std::istream &in, bool is_an_axiom, bool use_metric);
    ExplicitOperator(int cost, const std::string &name, bool is_an_axiom);

    // We ignore costs!
    bool operator==(const ExplicitOperator &other) const;
    bool operator<(const ExplicitOperator &other) const;

    int get_encoding_size() const;
};

class RootTask : public AbstractTask {
    friend class DecoupledRootTask;
protected:
    std::vector<ExplicitVariable> variables;
    // TODO: think about using hash sets here.
    std::vector<std::vector<std::set<FactPair>>> mutexes;
    std::vector<ExplicitOperator> operators;
    std::vector<ExplicitOperator> axioms;
    std::vector<int> initial_state_values;
    std::vector<FactPair> goals;

    const ExplicitVariable &get_variable(int var) const;
    virtual const ExplicitEffect &get_effect(int op_id, int effect_id, bool is_axiom) const;
    const ExplicitOperator &get_operator_or_axiom(int index, bool is_axiom) const;

public:
    explicit RootTask(std::istream &in);
    RootTask() {}

    virtual int get_num_variables() const override;
    virtual std::string get_variable_name(int var) const override;
    virtual int get_variable_domain_size(int var) const override;
    virtual int get_variable_axiom_layer(int var) const override;
    virtual int get_variable_default_axiom_value(int var) const override;
    virtual std::string get_fact_name(const FactPair &fact) const override;
    virtual bool are_facts_mutex(
        const FactPair &fact1, const FactPair &fact2) const override;

    virtual int get_operator_cost(int index, bool is_axiom) const override;
    virtual std::string get_operator_name(
        int index, bool is_axiom) const override;
    virtual int get_num_operators() const override;
    virtual int get_num_operator_preconditions(
        int index, bool is_axiom) const override;
    virtual FactPair get_operator_precondition(
        int op_index, int fact_index, bool is_axiom) const override;
    virtual int get_num_operator_effects(
        int op_index, bool is_axiom) const override;
    virtual int get_num_operator_effect_conditions(
        int op_index, int eff_index, bool is_axiom) const override;
    virtual FactPair get_operator_effect_condition(
        int op_index, int eff_index, int cond_index, bool is_axiom) const override;
    virtual FactPair get_operator_effect(
        int op_index, int eff_index, bool is_axiom) const override;
    virtual int convert_operator_index(
        int index, const AbstractTask *ancestor_task) const override;

    virtual int get_num_axioms() const override;

    virtual int get_num_goals() const override;
    virtual FactPair get_goal_fact(int index) const override;

    virtual std::vector<int> get_initial_state_values() const override;
    virtual void convert_ancestor_state_values(
        std::vector<int> &values,
        const AbstractTask *ancestor_task) const override;

    virtual int get_encoding_size(bool with_mutexes) const;
};
}
#endif

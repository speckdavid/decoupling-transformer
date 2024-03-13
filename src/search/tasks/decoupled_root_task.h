#ifndef TASKS_DECOUPLED_ROOT_TASK_H
#define TASKS_DECOUPLED_ROOT_TASK_H

#include "root_task.h"

#include "../task_proxy.h"

#include <map>
#include <unordered_map>

namespace plugins {
class Options;
}

namespace decoupling {
class Factoring;
}

namespace tasks {
enum ConclusiveLeafEncoding {BASIC = 0, BINARY = 1, MULTIVALUED = 2};

/*
  Task transformation that decoupled the search space using derived variables and axioms
*/
class DecoupledRootTask : public RootTask {
    std::shared_ptr<RootTask> original_root_task;
    std::shared_ptr<decoupling::Factoring> factoring;

    bool skip_unnecessary_leaf_effects;
    bool same_leaf_preconditons_single_variable;
    bool conclusive_operators;
    ConclusiveLeafEncoding conclusive_leaf_encoding;
    mutable std::unordered_map<int, std::unordered_set<int>> original_operator_tr_eff_vars;

    std::unordered_map<int, int> center_var_to_pvar;
    std::unordered_map<int, int> conclusive_leaf_var_to_pvar;
    std::unordered_map<int, int> leaf_to_goal_svar;
    std::unordered_map<int, std::unordered_map<int, int>> leaf_lstate_to_pvar;
    std::unordered_map<int, std::unordered_map<int, int>> leaf_lstate_to_svar;

    std::map<std::vector<FactPair>, int> precondition_to_svar;
    std::unordered_map<int, std::unordered_map<int, int>> leaf_op_to_svar;

    std::unordered_map<int, int> global_op_id_to_original_op_id;
    std::unordered_map<int, int> original_op_id_to_global_op_id;

public:
    DecoupledRootTask(const plugins::Options &options);
    virtual ~DecoupledRootTask() override = default;

    int get_original_operator_id(int op_id) const;

    void reconstruct_plan_if_necessary(std::vector<OperatorID> &path,
                                       std::vector<State> &states) const override;

    virtual TaskProxy get_task_proxy_for_plan_saving() const {
        // If we run decoupled search, we need the original task to save the reconstructed plan.
        return TaskProxy(*original_root_task);
    }

    bool is_valid_decoupled_state(const State &dec_state) const;

    std::shared_ptr<AbstractTask> get_original_root_task() const;

protected:
    void print_statistics() const;
    void write_sas_file(const std::string &file_name) const;
    void write_pddl_files(const std::string &domain_file_name, const std::string &problem_file_name) const;
    void write_factoring_file(const std::string &file_name) const;

    bool are_initial_states_consistent() const;
    bool is_conclusive_operator(int op_id, int leaf) const;
    bool is_conclusive_leaf(int leaf) const;

    // variables
    std::vector<std::string> get_fact_names(const std::string &var_name) const;
    void create_center_variables();
    void create_leaf_state_variables();
    void create_goal_condition_variables();
    void create_precondition_variables();
    void create_variables();

    void create_mutexes();
    void create_initial_state();
    void create_goal();

    // operators
    void set_preconditions_of_operator(int op_id, ExplicitOperator &op);
    void set_center_effects_of_operator(int op_id, ExplicitOperator &op);
    void set_general_leaf_effects_of_operator(int op_id, ExplicitOperator &op, int leaf);
    void set_conclusive_leaf_effects_of_operator(int op_id, ExplicitOperator &op, int leaf, ConclusiveLeafEncoding encoding);
    void set_leaf_effects_of_operator(int op_id, ExplicitOperator &op);
    void create_operator(int op_id);
    void create_operators();

    // axioms
    void create_frame_axioms();
    void create_goal_axioms();
    void create_precondition_axioms();
    void create_leaf_only_operator_axioms();
    void create_axioms();

    void release_memory();

    void dump() const;
};
}

#endif

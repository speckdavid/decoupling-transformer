#ifndef DECOUPLING_FACTORING
#define DECOUPLING_FACTORING

#include "leaf_state_id.h"
#include "leaf_state_space.h"
#include "interaction_graph.h"

#include "../abstract_task.h"
#include "../task_proxy.h"

#include "../plugins/options.h"
#include "../utils/countdown_timer.h"

#include <memory>
#include <set>

namespace decoupling {

class InteractionGraph;

class Factoring : public std::enable_shared_from_this<Factoring> {
    friend class LeafStateSpace;
    friend class PathPrices;

    bool is_factoring_possible() const;
    bool is_two_leaf_factoring_possible() const;
    void check_factoring() const;

    bool optimize_leaf_unique_lstate;
    bool prune_fork_leaf_state_spaces;

    std::unique_ptr<LeafStateSpace> leaf_state_space;

    std::unique_ptr<InteractionGraph> interaction_graph;

    std::vector<std::vector<FactPair>> goals_by_leaf;

    std::vector<FactorID> var_to_factor;
    std::vector<int> var_to_id_in_factor;

    int num_global_operators;
    std::vector<bool> is_global_operator_;
    std::vector<std::vector<bool>> has_op_leaf_pre;
    std::vector<std::vector<bool>> has_op_leaf_eff;
    std::vector<std::vector<OperatorID>> leaf_operators;

    std::vector<bool> can_optimize_leaf_unique_lstate;

    void remove_never_applicable_global_ops(FactorID leaf);

    bool does_op_uniquely_fix_lstate(OperatorProxy op, FactorID leaf) const;
    bool does_op_restrict_leaf(OperatorProxy op, FactorID leaf) const;

    void check_can_optimize_leaf_unique_lstate();

    const std::vector<OperatorID> &get_leaf_operators(FactorID leaf) const;

    bool has_leaf_goal(FactorID leaf) const;

    bool is_center_applicable(const State &state, OperatorProxy op) const;
    int get_num_effects_on_leaf(OperatorProxy op, FactorID leaf) const;

    const std::vector<FactPair> &get_leaf_goals(FactorID factor) const;

protected:
    mutable utils::LogProxy log;
    utils::CountdownTimer factoring_timer;

    std::shared_ptr<AbstractTask> task;
    TaskProxy task_proxy;

    std::vector<int> center;
    std::vector<std::vector<int>> leaves;

    struct ActionSchema {
        int num_actions; // number of actions with the action schema
        std::vector<int> pre_vars; // sorted
        std::vector<int> eff_vars; // sorted

        ActionSchema(int num_actions, const std::vector<int> &pre_vars, const std::vector<int> &eff_vars)
            : num_actions(num_actions), pre_vars(pre_vars), eff_vars(eff_vars) {
        }

        void inc_num_actions() {
            num_actions++;
        }
    };

    std::vector<ActionSchema> action_schemas;

    std::vector<std::set<int>> var_to_affecting_op;

    int min_number_leaves;
    int max_leaf_size;

    void compute_action_schemas();
    void compute_var_to_ops_map();
    bool check_timeout() const;

    explicit Factoring(const plugins::Options &opts);
    void apply_factoring();
    void print_factoring() const;

    virtual void compute_factoring_() = 0;
    virtual void save_memory();

public:
    virtual ~Factoring() = default;

    void compute_factoring();

    bool is_center_variable(int var) const;
    bool is_leaf_variable(int var) const;
    int get_factor(int var) const;
    int get_id_in_factor(int var) const;

    bool is_leaf_only_operator(int operator_id) const;
    bool is_global_operator(int operator_id) const;

    int get_num_leaves() const;
    int get_num_leaf_variables(int leaf) const;
    int get_num_leaf_states(int leaf) const;
    int get_num_all_leaf_states() const;
    int get_num_all_goal_leaf_states() const;
    int get_num_global_operators() const;

    bool is_fork_leaf(FactorID leaf) const;
    bool is_ifork_leaf(FactorID leaf) const;
    bool is_fork_factoring() const;

    bool has_pre_on_leaf(OperatorID op_id, FactorID leaf) const;
    bool has_pre_on_leaf(int, int leaf) const;
    bool has_eff_on_leaf(OperatorID op_id, FactorID leaf) const;
    bool has_eff_on_leaf(int op_id, int leaf) const;
    bool has_pre_or_eff_on_leaf(OperatorID op_id, FactorID leaf) const;
    bool has_pre_or_eff_on_leaf(int op_id, int leaf) const;

    const std::vector<int> &get_center() const;
    const std::vector<std::vector<int>> &get_leaves() const;
    const std::vector<int> &get_leaf(int leaf) const;

    int get_initial_leaf_state(int leaf) const;
    const std::vector<LeafStateHash> &get_goal_leaf_states(int leaf) const;

    bool is_conclusive_leaf(FactorID leaf) const;
    bool is_conclusive_leaf(int leaf) const;

    std::vector<FactPair> get_leaf_state_values(int leaf, int leaf_state) const;

    // NOTE: this function is not very efficiently implemented and should be called sparsely
    std::vector<int> get_valid_leaf_states(int leaf, const std::vector<FactPair>& partial_state);

    std::vector<int> get_predecessors(int leaf, int leaf_state, int operator_id) const;

    void add_leaf_facts_to_state(std::vector<int> &state, int leaf, int leaf_state) const;

    std::string get_leaf_name(int leaf) const;
    std::string get_leaf_state_name(int leaf, int leaf_state) const;

    void insert_leaf_paths(std::vector<OperatorID> &path,
                           std::vector<State> &states,
                           const std::shared_ptr<AbstractTask> &original_root_task) const;

    static void add_options_to_feature(plugins::Feature &feature);
};
}

#endif

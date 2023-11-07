#ifndef DECOUPLING_FACTORING
#define DECOUPLING_FACTORING

#include "../abstract_task.h"
#include "../task_proxy.h"

#include "../plugins/options.h"
#include "../utils/countdown_timer.h"

#include <memory>
#include <set>

namespace decoupling {

class Factoring {
    utils::CountdownTimer factoring_timer;

protected:

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

        void incr_num_action() {
            num_actions++;
        }
    };

    static std::vector<ActionSchema> action_schemas;

    static std::vector<std::set<int> > var_to_affecting_op;

    int min_number_leaves;
    int max_leaf_size;

    void compute_action_schemas();
    void compute_var_to_ops_map();
    bool check_timeout() const;

    explicit Factoring(const plugins::Options &opts);
    void apply_factoring();
    void print_factoring() const;

    // TODO find a better name for this
    virtual void compute_factoring_() = 0;

public:
    virtual ~Factoring() = default;

    void compute_factoring();

    static void add_options_to_feature(plugins::Feature &feature);
};
}

#endif

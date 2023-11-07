#ifndef DECOUPLING_LP_FACTORING
#define DECOUPLING_LP_FACTORING

#include "factoring.h"

#include "../lp/lp_solver.h"

#include "../plugins/options.h"

#include <set>

namespace decoupling {

enum class STRATEGY {
    MML, // maximize mobile leaves
    MMAS, // maximize mobile action schemas
    MM_OPT, // maximize mobility
    MM_APPROX,
    MFA, // maximize mobile facts
    MM // maximize mobility (sum)
};

class LPFactoring : public decoupling::Factoring {

    struct PotentialLeaf {
        size_t id;
        int num_actions; // number of actions with the effect schema
        double max_flexibility;
        std::vector<int> vars; // sorted
        std::vector<size_t> action_schemes;
        std::vector<size_t> self_mobile_as; // subset of leaf_only_scheme whose pre_vars are in vars
        std::vector<double> as_flexibility;

        PotentialLeaf(size_t id, const std::vector<int> &vars)
        : id(id), num_actions(0), max_flexibility(0.0), vars(vars) {}

        void add_leaf_only_schema(int action_schema);
    };

    double infty;

    STRATEGY strategy;

    int min_mobility;

    double min_flexibility;

    double min_fact_flexibility;

    bool add_cg_sccs_;

    int max_merge_steps;

    bool merge_overlapping;

    bool merge_dependent;

    bool ignore_center_preconditions;

    std::vector<std::set<int>> var_to_p_leaves; // maps variables to potential leaf ids

    std::vector<PotentialLeaf> potential_leaves;

    void filter_potential_leaves();

    void compute_leaf_flexibility();

    std::vector<PotentialLeaf> compute_potential_leaves();

    void add_cg_sccs();

    void merge_potential_leaves();

    void recompute_var_to_p_leaves();

    void compute_factoring();

public:

    explicit LPFactoring(const plugins::Options &opts);

    static void add_options_to_parser(plugins::Feature &feature);

};
}

#endif

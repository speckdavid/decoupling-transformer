#ifndef DECOUPLING_LP_FACTORING
#define DECOUPLING_LP_FACTORING

#include "factoring.h"

#include "../algorithms/named_vector.h"

#include <set>

namespace lp {
class LPVariable;
class LPConstraint;
}

namespace plugins {
class Options;
class Feature;
}

namespace decoupling {

enum class STRATEGY {
    MML, // maximize mobile leaves
    MMAS, // maximize mobile action schemas
    MM_OPT, // maximize mobility
    MM_APPROX,
    MFA, // maximize mobile facts
    MM, // maximize mobility (sum)
    MCL, // maximize number of mobile conclusive leaves
    MCM, // maximize conclusive mobility, i.e. number of conclusive actions
};

class LPFactoring : public decoupling::Factoring {

    struct PotentialLeaf {
        const LPFactoring *factoring;
        size_t id;
        int num_actions; // number of actions with the effect schema
        double max_flexibility;
        std::vector<int> vars; // sorted
        std::vector<size_t> action_schemes;
        std::vector<size_t> self_mobile_as; // subset of action_schemes whose pre_vars are in vars
        std::vector<double> as_flexibility;

        PotentialLeaf(const LPFactoring *factoring, size_t id, const std::vector<int> &vars)
        : factoring(factoring), id(id), num_actions(0), max_flexibility(0.0), vars(vars) {}

        void add_leaf_only_schema(size_t action_schema);
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

    std::vector<std::set<int>> var_to_p_leaves; // maps variables to potential leaf ids

    std::vector<PotentialLeaf> potential_leaves;

    void filter_potential_leaves();

    void compute_leaf_flexibility();

    void compute_potential_leaves();

    void add_cg_sccs();

    void merge_potential_leaves();

    void recompute_var_to_p_leaves();

    std::vector<size_t> add_center_variables_and_get_ids(named_vector::NamedVector<lp::LPVariable> &variables,
                                                         named_vector::NamedVector<lp::LPConstraint> &constraints,
                                                         const std::vector<bool> &can_be_leaf_var) const;

    void add_leaf_intersection_constraints(named_vector::NamedVector<lp::LPConstraint> &constraints) const;

    void add_min_flexibility_constraints(named_vector::NamedVector<lp::LPConstraint> &constraints,
                                         const std::vector<std::vector<size_t>> &mob_as_var_ids) const;

    void add_min_mobility_constraints(named_vector::NamedVector<lp::LPConstraint> &constraints,
                                      const std::vector<std::vector<size_t>> &mob_as_var_ids) const;

    void add_potential_leaf_to_action_schema_constraints(named_vector::NamedVector<lp::LPConstraint> &constraints,
                                                         const std::vector<std::vector<size_t>> &mob_as_var_ids,
                                                         const std::vector<bool> &can_be_leaf_var,
                                                         const std::vector<size_t> &c_vars_ids) const;

    bool is_as_leaf_irrelevant(const ActionSchema &as, const PotentialLeaf &leaf) const;

    bool is_as_leaf_conclusive(const ActionSchema &as, const PotentialLeaf &leaf) const;

    bool has_as_pre_or_eff_on_leaf(const ActionSchema &as, const PotentialLeaf &leaf) const;

    void construct_lp_conclusive(named_vector::NamedVector<lp::LPVariable> &variables,
                                 named_vector::NamedVector<lp::LPConstraint> &constraints);

    void construct_lp_all(named_vector::NamedVector<lp::LPVariable> &variables,
                          named_vector::NamedVector<lp::LPConstraint> &constraints);

    void compute_factoring_() override;

    void save_memory() override;

public:

    explicit LPFactoring(const plugins::Options &opts);

    static void add_options_to_parser(plugins::Feature &feature);

};
}

#endif

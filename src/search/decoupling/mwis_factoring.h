#ifndef DECOUPLING_MWIS_FACTORING_H
#define DECOUPLING_MWIS_FACTORING_H

#include "factoring.h"


namespace plugins {
class Options;
class Feature;
}

namespace decoupling {

enum class WMIS_STRATEGY {
    MML, // maximize mobile leaves
    MMAS, // maximize mobile action schemas
    MM_OPT, // maximize mobility
    MFA, // maximize mobile facts
    MM, // maximize mobility (sum)
    MCL, // maximize number of mobile conclusive leaves
    MCM, // maximize conclusive mobility, i.e. number of conclusive actions
};

class MWISFactoring : public decoupling::Factoring {

    typedef std::vector<std::vector<int>> Graph;

    struct PotentialLeafNode {
        std::vector<int> outside_pre_vars;
        std::vector<int> vars;
        double weight;

        void dump(TaskProxy task_proxy) const {
            std::cout << " --------" << std::endl;
            std::cout << "outside pre vars:" << std::endl;
            for (int var : outside_pre_vars){
                std::cout << "\t" << task_proxy.get_variables()[var].get_fact(0).get_name() << std::endl;
            }
            std::cout << "vars:" << std::endl;
            for (int var : vars){
                std::cout << "\t"  << task_proxy.get_variables()[var].get_fact(0).get_name() << std::endl;
            }
            std::cout << "-------------------------" << std::endl;
        }

        bool operator==(const PotentialLeafNode &other) const {
            return vars == other.vars && outside_pre_vars == other.outside_pre_vars;
        }
    };

    struct PotentialLeaf {
        mutable int num_affecting_actions; // number of action have any effect on vars
        int num_actions; // number of actions with the effect schema
        std::vector<int> vars; // sorted
        std::vector<size_t> action_schemes;
        std::vector<size_t> self_mobile_as; // subset of action_schemes whose pre_vars are in vars

        explicit PotentialLeaf(const std::vector<int> &vars)
                : num_affecting_actions(0), num_actions(0), vars(vars) {
            assert(std::is_sorted(vars.begin(), vars.end()));
        }

        void add_leaf_only_schema(size_t as_id, const ActionSchema &action_schema);
    };

    WMIS_STRATEGY strategy;

    int min_mobility;

    double min_flexibility;

    double min_fact_flexibility;

    bool add_cg_sccs_;

    std::vector<std::vector<size_t>> variables_to_action_schemas;

    std::vector<PotentialLeafNode> potential_leaf_nodes; // final leaf candidates

    void compute_variables_to_action_schemas_map();

    static std::vector<size_t> get_superset_pleaf_ids(const PotentialLeaf &pleaf,
                                                      const std::vector<PotentialLeaf> &potential_leaves,
                                                      const std::vector<std::vector<size_t>> &var_to_p_leaves);

    void compute_fact_flexibility(
            std::vector<std::vector<std::unordered_map<size_t, int>>> &facts_to_mobility,
            std::vector<std::vector<int>> &sum_fact_mobility);

    void compute_potential_leaves();

    void add_cg_sccs(std::vector<PotentialLeaf> &potential_leaves,
                     std::vector<std::vector<size_t>> &var_to_p_leaves);

    void add_leaf_intersection_edges(Graph &graph,
                                     const std::vector<std::vector<size_t>> &var_to_p_leaves) const;

    void add_outside_pre_var_edges(Graph &graph,
                                   const std::vector<std::vector<size_t>> &var_to_p_leaves) const;

    bool fulfills_min_flexibility_and_mobility(
            const PotentialLeaf &pleaf,
            const std::vector<size_t> &included_as,
            const std::vector<std::vector<std::unordered_map<size_t, int>>> &facts_to_mobility,
            const std::vector<std::vector<int>> &sum_fact_mobility) const;

    void multiply_out_potential_leaf(const std::vector<std::pair<std::vector<int>, std::vector<size_t>>> &outside_pre_and_ases,
                                     const PotentialLeaf &pleaf,
                                     std::vector<int> &outside_pre_vars,
                                     std::vector<size_t> &included_as,
                                     size_t depth,
                                     int &ignored_leaf_candidates,
                                     const std::vector<std::vector<std::unordered_map<size_t, int>>> &facts_to_mobility,
                                     const std::vector<std::vector<int>> &sum_fact_mobility);

    void multiply_out_action_schemas(const std::vector<PotentialLeaf> &potential_leaves,
                                     const std::vector<std::vector<std::unordered_map<size_t, int>>> &facts_to_mobility,
                                     const std::vector<std::vector<int>> &sum_fact_mobility);

    bool is_as_leaf_irrelevant(const ActionSchema &as, const PotentialLeaf &leaf) const;

    static bool is_as_leaf_conclusive(const ActionSchema &as, const PotentialLeaf &leaf) ;

    static bool has_as_pre_or_eff_on_leaf(const ActionSchema &as, const PotentialLeaf &leaf) ;

    void construct_graph_conclusive_leaves(Graph &graph);

    void construct_graph(Graph &graph);

    std::vector<int> solve_wmis(const Graph &graph,
                                const std::vector<double> &weights,
                                const utils::CountdownTimer &timer);

    void compute_factoring_() override;

    void save_memory() override;

public:

    explicit MWISFactoring(const plugins::Options &opts);

    static void add_options_to_parser(plugins::Feature &feature);

};
}

#endif

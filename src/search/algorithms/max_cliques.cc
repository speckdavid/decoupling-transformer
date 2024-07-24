#include "max_cliques.h"

#include "../algorithms/dynamic_bitset.h"
#include "../utils/collections.h"
#include "../utils/countdown_timer.h"
#include "../utils/logging.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

using namespace std;

namespace max_cliques {
static vector<vector<int>> create_complement(const vector<vector<int>> &graph) {
    int n = graph.size();
    vector<vector<int>> complement_graph(n);

    dynamic_bitset::DynamicBitset is_neighbour(n);

    // Create the complement graph
    for (int i = 0; i < n; ++i) {
        is_neighbour.reset();
        for (int neighb : graph[i]){
            is_neighbour.set(neighb);
        }
        for (int j = i + 1; j < n; ++j) {
            if (!is_neighbour[j]) {
                complement_graph[i].push_back(j);
                complement_graph[j].push_back(i);
            }
        }
    }

    return complement_graph;
}

class MaxCliqueComputer {
    const vector<vector<int>> &graph;
    vector<vector<int>> &max_cliques;
    vector<int> current_max_clique;

    int get_maximizing_vertex(
        const vector<int> &subg, const vector<int> &cand) {
        assert(utils::is_sorted_unique(subg));
        assert(utils::is_sorted_unique(cand));

        //utils::g_log << "subg: " << subg << endl;
        //utils::g_log << "cand: " << cand << endl;
        size_t max = 0;
        // We will take the first vertex if there is no better one.
        int vertex = subg[0];

        for (size_t i = 0; i < subg.size(); ++i) {
            vector<int> intersection;
            intersection.reserve(subg.size());
            // for vertex u in subg get u's adjacent vertices: graph[subg[i]];
            set_intersection(cand.begin(), cand.end(),
                             graph[subg[i]].begin(), graph[subg[i]].end(),
                             back_inserter(intersection));

            if (intersection.size() > max) {
                max = intersection.size();
                vertex = subg[i];
                //utils::g_log << "success: there is a maximizing vertex." << endl;
            }
        }
        return vertex;
    }

    void expand(vector<int> &subg, vector<int> &cand) {
        // utils::g_log << "subg: " << subg << endl;
        // utils::g_log << "cand: " << cand << endl;
        if (subg.empty()) {
            //utils::g_log << "clique" << endl;
            max_cliques.push_back(current_max_clique);
        } else {
            int u = get_maximizing_vertex(subg, cand);

            vector<int> ext_u;
            ext_u.reserve(cand.size());
            set_difference(cand.begin(), cand.end(),
                           graph[u].begin(), graph[u].end(),
                           back_inserter(ext_u));

            while (!ext_u.empty()) {
                int q = ext_u.back();
                ext_u.pop_back();
                //utils::g_log << q << ",";
                current_max_clique.push_back(q);

                // subg_q = subg n gamma(q)
                vector<int> subg_q;
                subg_q.reserve(subg.size());
                set_intersection(subg.begin(), subg.end(),
                                 graph[q].begin(), graph[q].end(),
                                 back_inserter(subg_q));

                // cand_q = cand n gamma(q)
                vector<int> cand_q;
                cand_q.reserve(cand.size());
                set_intersection(cand.begin(), cand.end(),
                                 graph[q].begin(), graph[q].end(),
                                 back_inserter(cand_q));
                expand(subg_q, cand_q);

                // remove q from cand --> cand = cand - q
                cand.erase(lower_bound(cand.begin(), cand.end(), q));

                //utils::g_log << "back" << endl;
                current_max_clique.pop_back();
            }
        }
    }

public:
    MaxCliqueComputer(const vector<vector<int>> &graph_,
                      vector<vector<int>> &max_cliques_)
        : graph(graph_), max_cliques(max_cliques_) {
    }

    void compute() {
        vector<int> vertices_1;
        vertices_1.reserve(graph.size());
        for (size_t i = 0; i < graph.size(); ++i) {
            vertices_1.push_back(i);
        }
        vector<int> vertices_2(vertices_1);
        current_max_clique.reserve(graph.size());
        expand(vertices_1, vertices_2);
    }
};


void compute_max_cliques(
    const vector<vector<int>> &graph,
    vector<vector<int>> &max_cliques) {
    MaxCliqueComputer clique_computer(graph, max_cliques);
    clique_computer.compute();
}

void compute_max_independent_sets(
    const vector<vector<int>> &graph,
    vector<vector<int>> &max_independent_sets) {
    auto complement_graph = create_complement(graph);
    MaxCliqueComputer clique_computer(complement_graph, max_independent_sets);
    clique_computer.compute();
}

class MaxWeightCliqueComputer {
private:
    const vector<vector<int>> &graph;
    const vector<double> &node_weights;
    utils::CountdownTimer timer;

    // Result
    vector<int> incumbent_nodes;
    double incumbent_weight;

    void update_incumbent_if_improved(const vector<int> &C, double C_weight) {
        assert(utils::all_values_unique(C));
        if (C_weight > incumbent_weight) {
            incumbent_nodes = C;
            incumbent_weight = C_weight;

            #ifndef NDEBUG
            sort(incumbent_nodes.begin(), incumbent_nodes.end());
            utils::g_log << "New weight: " << incumbent_weight << " for "
                         << incumbent_nodes << endl;
            #endif
        }
    }

    vector<int> greedily_find_independent_set(vector<int> P) {
        vector<int> independent_set;
        while (!P.empty()) {
            int v = P[0];
            independent_set.push_back(v);
            P.erase(P.begin());
            P.erase(remove_if(P.begin(), P.end(), [this, v](int w) {
                                  return find(graph.at(v).begin(), graph.at(v).end(), w) != graph.at(v).end();
                              }), P.end());
        }
        return independent_set;
    }

    vector<int> find_branching_nodes(vector<int> P, double target) {
        // TODO: this method seems to be a bottleneck; probably worth switching to faster data structure
        unordered_map<int, double> residual_wt;
        for (int v : P) {
            residual_wt[v] = node_weights[v];
        }

        double total_wt = 0;
        while (!P.empty()) {
            vector<int> independent_set = greedily_find_independent_set(P);
            int min_element_in_class = *min_element(independent_set.begin(), independent_set.end(),
                                                    [&residual_wt](int a, int b) {
                                                        return residual_wt[a] < residual_wt[b];
                                                    });
            double min_wt_in_class = residual_wt[min_element_in_class];

            total_wt += min_wt_in_class;
            if (total_wt > target)
                break;

            for (int v : independent_set) {
                residual_wt[v] -= min_wt_in_class;
            }

            P.erase(remove_if(P.begin(), P.end(), [&residual_wt](int v) {
                                  return residual_wt[v] == 0;
                              }), P.end());
        }
        return P;
    }

    void expand(const vector<int> &C, double C_weight, vector<int> P) {
        update_incumbent_if_improved(C, C_weight);
        vector<int> branching_nodes = find_branching_nodes(P, incumbent_weight - C_weight);
        assert(utils::all_values_unique(branching_nodes));

        while (!branching_nodes.empty()) {
            int v = branching_nodes.back();
            branching_nodes.pop_back();
            P.erase(remove(P.begin(), P.end(), v), P.end());

            vector<int> new_C = C;
            new_C.push_back(v);

            double new_C_weight = C_weight + node_weights[v];
            vector<int> new_P;
            for (int w : P) {
                if (find(graph.at(v).begin(), graph.at(v).end(), w) != graph.at(v).end()) {
                    new_P.push_back(w);
                }
            }
            if (timer.is_expired()) {
                return;
            }
            expand(new_C, new_C_weight, new_P);
        }
    }

public:
    MaxWeightCliqueComputer(const vector<vector<int>> &graph, const vector<double> &weights, const utils::CountdownTimer &max_time)
        : graph(graph), node_weights(weights), timer(max_time), incumbent_weight(0) {
        if (graph.size() != weights.size()) {
            cerr << "Number of nodes does not match number of weights!" << endl;
            utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
        }
        if (graph.empty()) {
            cerr << "Graph is empty!" << endl;
            utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
        }
    }

    double find_max_weight_clique(vector<int> &max_clique) {
        if (graph.empty()) {
            return 0;
        }

        vector<int> nodes(graph.size());
        iota(nodes.begin(), nodes.end(), 0);

        sort(nodes.begin(), nodes.end(), [this](int a, int b) {
                 return graph.at(a).size() > graph.at(b).size();
             });
        nodes.erase(remove_if(nodes.begin(), nodes.end(), [this](int v) {
                                  return node_weights[v] <= 0;
                              }), nodes.end());
        assert(utils::all_values_unique(nodes));
        expand({}, 0, nodes);
        max_clique = incumbent_nodes;
        sort(max_clique.begin(), max_clique.end());
        return incumbent_weight;
    }
};

double compute_max_weighted_clique(
    const vector<vector<int>> &graph,
    const vector<double> &weights,
    vector<int> &max_clique,
    double max_time) {
    MaxWeightCliqueComputer computer(graph, weights, utils::CountdownTimer(max_time));
    return computer.find_max_weight_clique(max_clique);
}

double compute_max_weighted_independent_set(
    const vector<vector<int>> &graph,
    const vector<double> &weights,
    vector<int> &independent_set,
    double max_time) {
    utils::CountdownTimer max_timer(max_time);
    auto complement_graph = create_complement(graph);
    MaxWeightCliqueComputer computer(complement_graph, weights, max_timer);
    return computer.find_max_weight_clique(independent_set);
}
}

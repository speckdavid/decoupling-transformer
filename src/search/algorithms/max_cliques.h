#ifndef ALGORITHMS_MAX_CLIQUES_H
#define ALGORITHMS_MAX_CLIQUES_H

#include <limits>
#include <vector>

namespace max_cliques {
/* Implementation of the Max Cliques algorithm by Tomita et al. See:

   Etsuji Tomita, Akira Tanaka and Haruhisa Takahashi, The Worst-Case
   Time Complexity for Generating All Maximal Cliques. Proceedings of
   the 10th Annual International Conference on Computing and
   Combinatorics (COCOON 2004), pp. 161-170, 2004.

   In the paper the authors use a compressed output of the cliques,
   such that the algorithm is in O(3^{n/3}).

   This implementation is in O(n 3^{n/3}) because the cliques are
   output explicitly. For a better runtime it could be useful to
   use bit vectors instead of vectors.
 */
extern void compute_max_cliques(
    const std::vector<std::vector<int>> &graph,
    std::vector<std::vector<int>> &max_cliques);

extern void compute_max_independent_sets(
    const std::vector<std::vector<int>> &graph,
    std::vector<std::vector<int>> &max_independent_sets);

/* Implementation of the Weighted Max Clique algorithm implemented
   in the python networkx package:

   The implementation is recursive, and therefore it may run into recursion
   depth issues if G contains a clique whose number of nodes is close to the
   recursion depth limit.

   At each search node, the algorithm greedily constructs a weighted independent
   set cover of part of the graph in order to find a small set of nodes on which
   to branch. The algorithm is very similar to the algorithm of Tavares et al. [1],
   other than the fact that the NetworkX version does not use bitsets. This style of
   algorithm for maximum weight clique (and maximum weight independent set, which is the
   same problem but on the complement graph) has a decades-long history. See Algorithm B of
   Warren and Hicks [2] and the references in that paper.

   [1] Tavares, W.A., Neto, M.B.C., Rodrigues, C.D., Michelon, P.: Um algoritmo de branch and
       bound para o problema da clique máxima ponderada. Proceedings of XLVII SBPO 1 (2015).
   [2] Warren, Jeffrey S, Hicks, Illya V.: Combinatorial Branch-and-Bound for the Maximum Weight
       Independent Set Problem. Technical Report, Texas A&M University (2016).

    We stop after max_time seconds.
 */
extern double compute_max_weighted_clique(
    const std::vector<std::vector<int>> &graph,
    const std::vector<double> &weights,
    std::vector<int> &max_clique,
    double max_time = std::numeric_limits<double>::infinity());

extern double compute_max_weighted_independent_set(
    const std::vector<std::vector<int>> &graph,
    const std::vector<double> &weights,
    std::vector<int> &independent_set,
    double max_time = std::numeric_limits<double>::infinity());
}
#endif

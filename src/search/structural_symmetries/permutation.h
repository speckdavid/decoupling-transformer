#ifndef STRUCTURAL_SYMMETRIES_PERMUTATION_H
#define STRUCTURAL_SYMMETRIES_PERMUTATION_H

#include <vector>

#include "group.h"

namespace tasks {
class SymmetricRootTask;
}

/*
  This class represents search symmetries, i.e., it only stores a mapping of
  variables and values, but not of operators of the planning task, since the
  mapping of these is not used for symmetry pruning techniques.
*/
namespace structural_symmetries {
class Permutation {
    friend class tasks::SymmetricRootTask;
public:
    explicit Permutation(const Group &group);

    Permutation(const Group &group, const unsigned int *full_perm);

    Permutation(const Group &group, const std::vector<int> &full_perm);

    Permutation(const Permutation &perm, bool invert = false);

    Permutation(const Permutation &perm1, const Permutation &perm2);

    ~Permutation() = default;

    bool identity() const;

    void print_cycle_notation() const;

    void print_affected_variables_by_cycles() const;

    int get_value(int ind) const {
        return value[ind];
    }

    void dump_var_vals() const;

    void dump() const;

    void dump_fdr() const;

    std::pair<int, int> get_new_var_val_by_old_var_val(const int var, const int val) const;

    bool replace_if_less(std::vector<int> &state) const;

    void replace_partial_state(std::vector<int> &state) const;

    bool replace_if_less_partial_state(std::vector<int> &state) const;

    const std::vector<int> &get_affected_vars() const {
        return vars_affected;
    }

    bool affects_variable(int var) const;

private:
    const Group &group;
    std::vector<int> value;
    std::vector<int> vars_affected;

    // Need to keep the connection between affected vars, ie which var goes into which.
    std::vector<int> from_vars;
    // Affected vars by cycles
    std::vector<std::vector<int> > affected_vars_cycles;

    void finalize();

    void _allocate();

    void _copy_value_from_permutation(const Permutation &perm);

    void _inverse_value_from_permutation(const Permutation &perm);

    bool is_less_partial_state(const std::vector<int> &state) const;
};
}
#endif

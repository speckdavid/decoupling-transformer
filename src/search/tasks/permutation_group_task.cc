#include "permutation_group_task.h"

#include "../structural_symmetries/group.h"
#include "../structural_symmetries/permutation.h"
#include "../utils/logging.h"
#include "../utils/system.h"
#include "../task_utils/task_properties.h"

using namespace std;
using namespace structural_symmetries;

namespace extra_tasks {
PermutationGroupTask::PermutationGroupTask(
    const shared_ptr<AbstractTask> &parent,
    const shared_ptr<Group> &group)
    : DelegatingTask(parent),
      group(group) {
    TaskProxy parent_task_proxy(*parent);
    task_properties::verify_no_axioms(parent_task_proxy);
    task_properties::verify_no_conditional_effects(parent_task_proxy);

    if (!group->is_initialized()) {
        group->compute_symmetries(TaskProxy(*parent));
        if (!group->has_symmetries()) {
            utils::g_log << "No symmetries found, aborting.." << endl;
            utils::exit_with(utils::ExitCode::SEARCH_UNSOLVED_INCOMPLETE);
        }
    }


    for (int g = 0; g < group->get_num_generators(); ++g) {
        vector<tasks::ExplicitEffect> effects;
        auto const &gen = group->get_permutation(g);
        for (int var : gen.get_affected_vars()) {
            for (int val = 0; val < get_variable_domain_size(var); ++val) {
                auto [new_var, new_val] = gen.get_new_var_val_by_old_var_val(var, val);
                effects.emplace_back(var, val, vector<FactPair>{FactPair(new_var, new_val)});
            }
        }
        generator_effects.push_back(effects);
    }
}

int PermutationGroupTask::get_operator_cost(int index, bool is_axiom) const {
    assert(!is_axiom);
    return 1;
}

string PermutationGroupTask::get_operator_name(int index, bool is_axiom) const {
    assert(!is_axiom);
    return "Generator " + to_string(index);
}

int PermutationGroupTask::get_num_operators() const {
    return group->get_num_generators();
}

int PermutationGroupTask::get_num_operator_preconditions(int index, bool is_axiom) const {
    assert(!is_axiom);
    return 0;
}

FactPair PermutationGroupTask::get_operator_precondition(int op_index, int fact_index, bool is_axiom) const {
    utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    return FactPair::no_fact;
}

int PermutationGroupTask::get_num_operator_effects(int op_index, bool is_axiom) const {
    assert(!is_axiom);
    return generator_effects.at(op_index).size();
}

int PermutationGroupTask::get_num_operator_effect_conditions(
    int op_index, int eff_index, bool is_axiom) const {
    assert(!is_axiom);
    return generator_effects.at(op_index).at(eff_index).conditions.size();
}

FactPair PermutationGroupTask::get_operator_effect_condition(
    int op_index, int eff_index, int cond_index, bool is_axiom) const {
    assert(!is_axiom);
    return generator_effects.at(op_index).at(eff_index).conditions.at(cond_index);
}

FactPair PermutationGroupTask::get_operator_effect(
    int op_index, int eff_index, bool is_axiom) const {
    assert(!is_axiom);
    return generator_effects.at(op_index).at(eff_index).fact;
}

int PermutationGroupTask::convert_operator_index_to_parent(int index) {
    return OperatorID::no_operator.get_index();
}
}

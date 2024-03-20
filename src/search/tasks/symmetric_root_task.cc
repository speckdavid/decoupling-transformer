#include "symmetric_root_task.h"

#include "../state_registry.h"

#include "../plugins/plugin.h"
#include "../structural_symmetries/group.h"
#include "../structural_symmetries/permutation.h"
#include "../task_utils/successor_generator.h"
#include "../task_utils/dump_sas_task.h"
#include "../task_utils/task_properties.h"
#include "../utils/rng.h"

#include <fstream>

using namespace std;
using namespace structural_symmetries;

namespace tasks {
SymmetricRootTask::SymmetricRootTask(const plugins::Options &options)
    : RootTask(),
      original_root_task(dynamic_pointer_cast<RootTask>(tasks::g_root_task)),
      group(options.get<shared_ptr<Group>>("symmetries")),
      empty_value_strategy(options.get<EmptyValueStrategy>("empty_value_strategy")),
      compute_perfect_canonical(options.get<bool>("compute_perfect_canonical")) {
    TaskProxy original_task_proxy(*original_root_task);
    task_properties::verify_no_axioms(original_task_proxy);
    task_properties::verify_no_conditional_effects(original_task_proxy);

    if (!group->is_initialized()) {
        group->compute_symmetries(original_task_proxy);
        if (!group->has_symmetries()){
            utils::g_log << "No symmetries found, aborting.." << endl;
            utils::exit_with(utils::ExitCode::SEARCH_UNSOLVED_INCOMPLETE);
        }
    }

    utils::Timer transformation_timer;

    // copy everything from original root task
    // TODO avoid the copies
    variables = original_root_task->variables;
    goals = original_root_task->goals;
    mutexes = original_root_task->mutexes;

    if (empty_value_strategy == EmptyValueStrategy::NONE) {
        base_state_for_op_permutation = vector<int>(get_num_variables(), -1);
    } else if (empty_value_strategy == EmptyValueStrategy::INIT){
        base_state_for_op_permutation = original_root_task->initial_state_values;
    } else if (empty_value_strategy == EmptyValueStrategy::RANDOM) {
        // TODO add RNG to options
        base_state_for_op_permutation = vector<int>(get_num_variables());
        utils::RandomNumberGenerator rng;
        for (int var = 0; var < get_num_variables(); ++var) {
            base_state_for_op_permutation[var] = rng.random(variables[var].domain_size);
        }
    } else {
        utils::g_log << "ERROR: unknown empty_value_strategy: " << empty_value_strategy << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }

    utils::g_log << "Starting symmetry transformation!" << endl;
    utils::g_log << "Creating new initial state..." << endl;
    create_initial_state();
    utils::g_log << "Creating new operators..." << endl;
    create_operators();

    if (options.get<bool>("dump_task")) {
        dump();
    }

    utils::g_log << "Time for symmetry transformation: " << transformation_timer << endl;

    if (options.get<bool>("dump_task")) {
        dump();
    }

    print_statistics();

    release_memory();

    if (options.get<bool>("write_sas_file")) {
        write_sas_file("dec_output.sas");
        utils::exit_with(utils::ExitCode::SEARCH_UNSOLVED_INCOMPLETE);
    }
}

void SymmetricRootTask::reconstruct_plan_if_necessary(vector<OperatorID> &path,
                                                      vector<State> &states,
                                                      StateRegistry &state_registry) const {

    TaskProxy original_task_proxy(*original_root_task);

    vector<RawPermutation> permutations;

    for (int i = 0; i < static_cast<int>(states.size()) - 1; ++i) {
        OperatorID op_id = path[i];
        State parent_state = states[i+1];
        State new_state = state_registry.get_successor_state(parent_state, original_task_proxy.get_operators()[op_id]);

        RawPermutation p;
        if (new_state.get_id() != states[i].get_id()){
            const ExplicitOperator &original_op = original_root_task->get_operator_or_axiom(op_id.get_index(), false);
            if (compute_perfect_canonical) {
                Permutation perm(group->get_perfect_canonical_permutation(get_state_for_operator_permutation(original_op)),
                                 true);
                p = perm.value;
            } else {
                Permutation perm(group->get_canonical_permutation(get_state_for_operator_permutation(original_op)),
                                 true);
                p = perm.value;
            }
        } else {
            p = group->new_identity_raw_permutation();
        }
        permutations.push_back(std::move(p));
    }
    // invert initial-state permutation to obtain original initial state
    Permutation rev_init_state_p(Permutation(*initial_state_permutation), true);
    permutations.push_back(rev_init_state_p.value);

    assert(states.size() == permutations.size());

    vector<RawPermutation> reverse_permutations;
    RawPermutation temp_p = group->new_identity_raw_permutation();
    while (permutations.begin() != permutations.end()) {
        const RawPermutation &p = permutations.back();
        temp_p = group->compose_permutations(p, temp_p);
        reverse_permutations.push_back(temp_p);
        permutations.pop_back();
    }
    for (size_t i = 0; i < states.size(); ++i){
        const RawPermutation &permutation = reverse_permutations[states.size() - i-1];
        states[i] = state_registry.permute_state(states[i],
                                                 Permutation(*group, permutation));
    }
    reverse_permutations.clear();
    path.clear();
    for (int i = states.size() - 1; i > 0; i--) {
        vector<OperatorID> applicable_ops;
        successor_generator::g_successor_generators[TaskProxy(*this)].generate_applicable_ops(states[i], applicable_ops);
        bool found = false;
        int min_cost_op=0;
        int min_cost=numeric_limits<int>::max();

        for (size_t o = 0; o < applicable_ops.size(); o++) {
            OperatorProxy op = original_task_proxy.get_operators()[applicable_ops[o]];
            State succ_state = state_registry.get_successor_state(states[i], op);
            if (succ_state.get_id() == states[i-1].get_id()) {
                found = true;
                if (op.get_cost() < min_cost) {
                    min_cost = op.get_cost();
                    min_cost_op = o;
                }
            }
        }
        if (!found) {
            utils::g_log << "No operator is found!!!" << endl
                         << "Cannot reach the state" << endl;
            task_properties::dump_pddl(states[i-1]);
            utils::g_log << endl << "From the state" << endl;
            task_properties::dump_pddl(states[i]);
            utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
        }
        path.push_back(applicable_ops[min_cost_op]);
    }
    reverse(path.begin(), path.end());
}

void SymmetricRootTask::print_statistics() const {
    utils::g_log << "Original task size: " << original_root_task->get_encoding_size(false) << endl;
    utils::g_log << "Original number of variables: " << original_root_task->variables.size() << endl;
    utils::g_log << "Original number of operators: " << original_root_task->get_num_operators() << endl;

    utils::g_log << "Task size: " << get_encoding_size(false) << endl;
    utils::g_log << "Number of variables: " << variables.size() << endl;
    utils::g_log << "Number of operators: " << get_num_operators() << endl;
}

void SymmetricRootTask::write_sas_file(const std::string &file_name) const {
    utils::Timer write_sas_file_timer;
    utils::g_log << "Writing to dec_output.sas..." << flush;
    std::ofstream output_file;
    output_file.open(file_name);
    dump_sas_task::dump_as_SAS(*this, output_file);
    utils::g_log << "done!" << endl;
    utils::g_log << "Time for writing sas file: " << write_sas_file_timer << endl;
}

void SymmetricRootTask::create_initial_state() {
    initial_state_values = original_root_task->initial_state_values;
    if (compute_perfect_canonical) {
        initial_state_permutation = make_unique<Permutation>(
                group->do_perfect_canonical_inplace_and_get_permutation(initial_state_values));
    } else {
        initial_state_permutation = make_unique<Permutation>(
                group->do_canonical_inplace_and_get_permutation(initial_state_values));
    }
}

inline void add_conditional_permuted_effects(tasks::ExplicitOperator &new_op,
                                             const Permutation &perm,
                                             int from_var,
                                             int domain_size) {
    for (int from_val = 0; from_val < domain_size; ++from_val){
        auto [to_var, to_val] = perm.get_new_var_val_by_old_var_val(from_var, from_val);
        assert(to_var == perm.get_new_var_val_by_old_var_val(from_var, 0).first);
        new_op.effects.push_back({to_var, to_val, {{from_var, from_val}}});
    }
}

inline void set_partial_state_from_action(vector<int> &state, const ExplicitOperator &op) {
    for (const auto &pre : op.preconditions){
        state[pre.var] = pre.value;
    }
    for (const auto &eff : op.effects){
        state[eff.fact.var] = eff.fact.value;
    }
}

vector<int> SymmetricRootTask::get_state_for_operator_permutation(const ExplicitOperator &op) const {
    vector<int> pre_eff_state(base_state_for_op_permutation);
    set_partial_state_from_action(pre_eff_state, op);
    return pre_eff_state;
}

void SymmetricRootTask::set_symmetry_effects_of_operator(int op_id, tasks::ExplicitOperator &new_op) {
    // TODO implement this
    //  try several variants:
    //  1) get canonical projected onto op's variables => done
    //  2) fill in remaining variables with initial-state values => done
    //  3) fill in remaining variables with random values (same for every action?) => done
    //  4) do a form of context splitting and introduce multiple copies of every action with different symmetries

    const auto &op = original_root_task->operators[op_id];

    vector<int> pre_eff_state(get_state_for_operator_permutation(op));

    unique_ptr<Permutation> perm;
    if (compute_perfect_canonical) {
        perm = make_unique<Permutation>(group->get_perfect_canonical_permutation(pre_eff_state));
    } else {
        perm = make_unique<Permutation>(group->get_canonical_permutation(pre_eff_state));
    }

    if (perm->identity()){
        // no symmetries
        new_op.effects = op.effects;
        return;
    }

    // check for which variables we added effects because of the permutation
    vector<bool> handled_var(original_root_task->get_num_variables(), false);
    pre_eff_state = vector<int>(variables.size(), -1);
    set_partial_state_from_action(pre_eff_state, op);

    auto &affected_vars_cycles = perm->affected_vars_cycles;
    for (size_t i = 0; i < affected_vars_cycles.size(); i++) {
        if (affected_vars_cycles[i].size() == 1) {
            int from_var = affected_vars_cycles[i][0];
            int from_val = pre_eff_state[from_var];
            if (from_val != -1){
                auto [to_var, to_val] = perm->get_new_var_val_by_old_var_val(from_var, from_val);
                // effect is fixed by original precondition or effect
                new_op.effects.push_back({to_var, to_val, {}});
                handled_var[to_var] = true;
                continue;
            }
            add_conditional_permuted_effects(new_op, *perm, from_var, variables[from_var].domain_size);
            auto [to_var, _] = perm->get_new_var_val_by_old_var_val(from_var, 0);
            handled_var[to_var] = true;
            continue;
        }
        // Remembering one value to be rewritten last
        int last_var = affected_vars_cycles[i][affected_vars_cycles[i].size() - 1];
        int last_val = pre_eff_state[last_var];

        for (int j = affected_vars_cycles[i].size() - 1; j > 0; j--) {
            // writing into variable affected_vars_cycles[i][j]
            int from_var = affected_vars_cycles[i][j - 1];
            int from_val = pre_eff_state[from_var];
            if (from_val != -1) {
                auto [to_var, to_val] = perm->get_new_var_val_by_old_var_val(from_var, from_val);
                assert(to_var == affected_vars_cycles[i][j]);
                // effect is fixed by original precondition or effect
                new_op.effects.push_back({to_var, to_val, {}});
                handled_var[to_var] = true;
                continue;
            }
            add_conditional_permuted_effects(new_op, *perm, from_var, variables[from_var].domain_size);
            auto [to_var, _] = perm->get_new_var_val_by_old_var_val(from_var, 0);
            handled_var[to_var] = true;
        }
        // writing the last one
        if (last_val != -1) {
            auto [to_var, to_val] = perm->get_new_var_val_by_old_var_val(last_var, last_val);
            // effect is fixed by original precondition or effect
            new_op.effects.push_back({to_var, to_val, {}});
            handled_var[to_var] = true;
            continue;
        }
        add_conditional_permuted_effects(new_op, *perm, last_var, variables[last_var].domain_size);
        auto [to_var, _] = perm->get_new_var_val_by_old_var_val(last_var, 0);
        handled_var[to_var] = true;
    }

    // copy original effects of variables not affected by the permutation
    for (const auto &eff: op.effects) {
        if (!handled_var[eff.fact.var]){
            new_op.effects.push_back(eff);
        }
    }
}

void SymmetricRootTask::create_operator(int op_id) {
    const auto &op = original_root_task->operators[op_id];

    ExplicitOperator new_op(op.cost, op.name, op.preconditions); // TODO avoid copying preconditions
    set_symmetry_effects_of_operator(op_id, new_op);
    assert(!new_op.effects.empty());

    operators.push_back(new_op);
}

void SymmetricRootTask::create_operators() {
    for (size_t op_id = 0; op_id < original_root_task->operators.size(); ++op_id) {
        create_operator(op_id);
    }
}

void SymmetricRootTask::release_memory() {
    // TODO: release memory
}

void SymmetricRootTask::dump() const {
    task_properties::dump_task(TaskProxy(*this), true, true);
}


shared_ptr<AbstractTask> SymmetricRootTask::get_original_root_task() const {
    return original_root_task;
}

class SymmetricRootTaskFeature : public plugins::TypedFeature<AbstractTask, SymmetricRootTask> {
public:
    SymmetricRootTaskFeature() : TypedFeature("symmetry") {
        document_title("Symmetric task");
        document_synopsis(
            "A symmetry transformation of the root task.");

        add_option<shared_ptr<Group>>("symmetries", "method that computes symmetries");
        add_option<EmptyValueStrategy>("empty_value_strategy", "How to treat variables not mentioned in operator precondition and effect", "none");
        add_option<bool>("compute_perfect_canonical", "Computes the perfect canonical for each orbit.", "true");
        add_option<bool>("write_sas_file", "Writes the decoupled task to dec_output.sas and terminates.", "false");
        add_option<bool>("dump_task", "Dumps the task to the console", "false");
    }

    virtual shared_ptr<SymmetricRootTask> create_component(const plugins::Options &options, const utils::Context &) const override {
        return make_shared<SymmetricRootTask>(options);
    }
};

static plugins::FeaturePlugin<SymmetricRootTaskFeature> _plugin;

static plugins::TypedEnumPlugin<EmptyValueStrategy> _enum_plugin({
        {"none", "skip variables not set in precondition or effect"},
        {"init", "use initial state values for variables not set in precondition or effect"},
        {"random", "use random values for variables not set in precondition or effect"}
    });
}

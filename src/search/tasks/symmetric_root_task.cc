#include "symmetric_root_task.h"

#include "permutation_group_task.h"

#include "../state_registry.h"

#include "../plugins/plugin.h"
#include "../structural_symmetries/group.h"
#include "../structural_symmetries/permutation.h"
#include "../task_utils/dump_sas_task.h"
#include "../task_utils/successor_generator.h"
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
          skip_mutex_preconditions(options.get<bool>("skip_mutex_preconditions")),
          skip_unaffected_variables(options.get<bool>("skip_unaffected_variables")),
          skip_unaffected_variables_relevant_permutations(options.get<bool>("skip_unaffected_variables_relevant_permutations")),
          context_splitting(options.get<bool>("context_splitting")),
          decoupled_splitting(options.get<bool>("decoupled_splitting")),
          max_number_contexts_per_operator(options.get<int>("max_number_contexts_per_operator")),
          compute_perfect_canonical(options.get<bool>("compute_perfect_canonical")) {
    TaskProxy original_task_proxy(*original_root_task);
    task_properties::verify_no_axioms(original_task_proxy);
    task_properties::verify_no_conditional_effects(original_task_proxy);

    if (!context_splitting){
        // TODO: print warning?
        skip_mutex_preconditions = false;
        skip_unaffected_variables = false;
        skip_unaffected_variables_relevant_permutations = false;
        decoupled_splitting = false;
        max_number_contexts_per_operator = 0;
    }

    if (!group->is_initialized()) {
        group->compute_symmetries(original_task_proxy);
        if (!group->has_symmetries()) {
            utils::g_log << "No symmetries found, aborting.." << endl;
            utils::exit_with(utils::ExitCode::SEARCH_UNSOLVED_INCOMPLETE);
        }
    }

    if (context_splitting && decoupled_splitting && group->get_permutation_components().size() == 1){
        // TODO evaluate which variant is actually better
        utils::g_log << "WARNING: permutation interaction graph is strongly connected, disabling decoupled_splitting." << endl;
        decoupled_splitting = false;
    }

    utils::Timer transformation_timer;

    // copy everything from original root task
    // TODO avoid the copies
    variables = original_root_task->variables;
    goals = original_root_task->goals;
    mutexes = original_root_task->mutexes;

    if (context_splitting && skip_unaffected_variables_relevant_permutations){
        compute_decoupled_splitting_implied_relevant_vars();
    }

    if (skip_mutex_preconditions){
        // check if there are any mutexes and disable redundant somewhat expensive checking if not
        skip_mutex_preconditions = false;
        for (const auto &var_m : mutexes){
            for (const auto &fact_m : var_m){
                if (!fact_m.empty()){
                    skip_mutex_preconditions = true;
                    break;
                }
            }
            if (skip_mutex_preconditions){
                break;
            }
        }
    }


    if (empty_value_strategy == NONE) {
        base_state_for_op_permutation = vector<int>(get_num_variables(), -1);
    } else if (empty_value_strategy == INIT) {
        base_state_for_op_permutation = original_root_task->initial_state_values;
    } else if (empty_value_strategy == RANDOM) {
        // TODO add RNG to options
        utils::RandomNumberGenerator rng;
        rng.seed(42);
        base_state_for_op_permutation = vector<int>(get_num_variables());
        for (int var = 0; var < get_num_variables(); ++var) {
            base_state_for_op_permutation[var] = rng.random(variables[var].domain_size);
        }
    } else if (empty_value_strategy == GOAL) {
        base_state_for_op_permutation = vector<int>(get_num_variables(), -1);
        for (const auto &goal : original_root_task->goals) {
            base_state_for_op_permutation[goal.var] = goal.value;
        }
    } else if (empty_value_strategy == INIT_GOAL) {
        base_state_for_op_permutation = original_root_task->initial_state_values;
        for (const auto &goal : original_root_task->goals) {
            base_state_for_op_permutation[goal.var] = goal.value;
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
        task_properties::dump_task(original_task_proxy, true, true);
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

void SymmetricRootTask::compute_decoupled_splitting_implied_relevant_vars() {
    // we store implied variables for one step here, not recursively,
    // since we don't know the post condition, which breaks the recursion.
    // that means that we still need to do a fixed-point computation for every operator
    // in get_split_variables(op)
    assert(decoupled_splitting_implied_relevant_vars.empty());
    decoupled_splitting_implied_relevant_vars.resize(RootTask::get_num_variables());
    for (int var = 0; var < RootTask::get_num_variables(); ++var) {
        for (const auto &perm: group->generators) {
            if (perm.affects_variable(var)) {
                for (int implied_var : perm.vars_affected){
                    if (var != implied_var) {
                        decoupled_splitting_implied_relevant_vars[var].push_back(implied_var);
                    }
                }
            }
        }
        utils::sort_unique(decoupled_splitting_implied_relevant_vars[var]);
    }
}

vector<int> SymmetricRootTask::get_split_variables(const ExplicitOperator &op) const {
    vector<int> split_vars;
    if (skip_unaffected_variables_relevant_permutations){
        // obtain all variables that are (possibly recursively) affected by a permutation
        // that affects a variable in the operator's post condition

        vector<bool> is_relevant_var(RootTask::get_num_variables(), false);
        for (const auto &pre : op.preconditions){
            is_relevant_var[pre.var] = true;
            for (int var : decoupled_splitting_implied_relevant_vars[pre.var]){
                is_relevant_var[var] = true;
            }
        }
        for (const auto &eff : op.effects){
            is_relevant_var[eff.fact.var] = true;
            for (int var : decoupled_splitting_implied_relevant_vars[eff.fact.var]){
                is_relevant_var[var] = true;
            }
        }

        vector<int> post_condition_state(get_operator_post_condition(op, false));

        // for full pruning power we need to recursively collect the permutations that affect var
        // and then go over all variables again to check for these new permutations
        bool change = true;
        while (change){
            change = false;
            vector<int> added_vars;
            for (int check_var = 0; check_var < RootTask::get_num_variables(); ++check_var) {
                if (!is_relevant_var[check_var]){
                    continue;
                }
                for (int var: decoupled_splitting_implied_relevant_vars[check_var]) {
                    if (post_condition_state[var] == -1 && !is_relevant_var[var]) {
                        is_relevant_var[var] = true;
                        added_vars.push_back(var);
                    }
                }
            }
            for (int added_var : added_vars){
                for (int var: decoupled_splitting_implied_relevant_vars[added_var]) {
                    if (post_condition_state[var] == -1 && !is_relevant_var[var]) {
                        is_relevant_var[var] = true;
                        change = true;
                    }
                }
            }
        }

        // split variables are affected by some relevant permutation but are not in post condition
        for (int var = 0; var < RootTask::get_num_variables(); ++var){
            if (is_relevant_var[var] && post_condition_state[var] == -1){
                split_vars.push_back(var);
            }
        }
    } else {
        // split variables are any variables that are affected by some permutation and are not in post condition
        vector<int> post_condition_state(get_operator_post_condition(op, false));
        for (int var = 0; var < RootTask::get_num_variables(); ++var){
            if (post_condition_state[var] == -1){
                if (skip_unaffected_variables) {
                    if (group->is_var_affected_by_permutation(var)) {
                        split_vars.push_back(var);
                    }
                } else {
                    split_vars.push_back(var);
                }
            }
        }
    }
    if (max_number_contexts_per_operator < numeric_limits<int>::max()){
        // TODO: try different variants:
        //  1) simply take first k variables according to FD variable order => done
        //  2) like 1) but inverse variable order
        //  3) try to maximize the number of completely covered symmetry components
        //  4) take some variables from all components
        // for partial splitting, sort split_vars here and truncate it accordingly
        utils::sort_unique(split_vars); // to match FD variable order
        int size = 1;
        for (size_t i = 0; i < split_vars.size(); ++i){
            size *= RootTask::get_variable_domain_size(split_vars[i]);
            if (size > max_number_contexts_per_operator){
                split_vars.resize(i);
                break;
            }
        }
    }
    return split_vars;
}

void SymmetricRootTask::reconstruct_plan_if_necessary(vector<OperatorID> &path,
                                                      vector<State> &states,
                                                      StateRegistry &state_registry) const {
    TaskProxy original_task_proxy(*original_root_task);

    vector<RawPermutation> permutations;

    for (int i = 0; i < static_cast<int>(states.size()) - 1; ++i) {
        OperatorID op_id = path[i];
        State &parent_state = states[i + 1];

        OperatorID original_op_id = op_id;
        if (!decoupled_splitting && context_splitting) {
            original_op_id = OperatorID(new_op_id_to_original_op_id[op_id.get_index()]);
        }

        State new_state = state_registry.get_successor_state(parent_state,
                                                             original_task_proxy.get_operators()[original_op_id]);

        RawPermutation p;
        if (new_state.get_id() != states[i].get_id()) {
            ExplicitOperator original_op = original_root_task->get_operator_or_axiom(original_op_id.get_index(), false);
            if (context_splitting) {
                if (decoupled_splitting){
                    vector<int> split_vars(get_split_variables(original_op));
                    for (int var : split_vars){
                        original_op.preconditions.emplace_back(var, parent_state[var].get_value());
                    }
                } else {
                    original_op.preconditions = operators[op_id.get_index()].preconditions;
                }
            }
            Permutation inv_perm(*get_permutation_for_operator(original_op), true);
            p = inv_perm.value;
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
    for (size_t i = 0; i < states.size(); ++i) {
        const RawPermutation &permutation = reverse_permutations[states.size() - i - 1];
        states[i] = state_registry.permute_state(states[i],
                                                 Permutation(*group, permutation));
    }
    reverse_permutations.clear();
    path.clear();
    for (int i = states.size() - 1; i > 0; i--) {
        vector<OperatorID> applicable_ops;
        successor_generator::g_successor_generators[TaskProxy(*this)].generate_applicable_ops(states[i],
                                                                                              applicable_ops);
        bool found = false;
        int min_cost_op = 0;
        int min_cost = numeric_limits<int>::max();

        for (size_t o = 0; o < applicable_ops.size(); o++) {
            OperatorID appl_op_id = applicable_ops[o];
            if (!decoupled_splitting && context_splitting) {
                appl_op_id = OperatorID(new_op_id_to_original_op_id[appl_op_id.get_index()]);
            }
            OperatorProxy op = original_task_proxy.get_operators()[appl_op_id];
            State succ_state = state_registry.get_successor_state(states[i], op);
            if (succ_state.get_id() == states[i - 1].get_id()) {
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
            task_properties::dump_pddl(states[i - 1]);
            utils::g_log << endl << "From the state" << endl;
            task_properties::dump_pddl(states[i]);
            utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
        }
        if (!decoupled_splitting && context_splitting) {
            path.push_back(OperatorID(new_op_id_to_original_op_id[applicable_ops[min_cost_op].get_index()]));
        } else {
            path.push_back(applicable_ops[min_cost_op]);
        }
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
    for (int from_val = 0; from_val < domain_size; ++from_val) {
        auto [to_var, to_val] = perm.get_new_var_val_by_old_var_val(from_var, from_val);
        assert(to_var == perm.get_new_var_val_by_old_var_val(from_var, 0).first);
        new_op.effects.push_back({to_var, to_val, {{from_var, from_val}}});
    }
}

inline void set_partial_state_from_action(vector<int> &state, const ExplicitOperator &op) {
    for (const auto &pre : op.preconditions) {
        state[pre.var] = pre.value;
    }
    for (const auto &eff : op.effects) {
        state[eff.fact.var] = eff.fact.value;
    }
}

vector<int> SymmetricRootTask::get_operator_post_condition(const ExplicitOperator &op,
                                                           bool fill_with_base_state) const {
    vector<int> pre_eff_state;
    if (fill_with_base_state) {
        pre_eff_state = base_state_for_op_permutation;
    } else {
        pre_eff_state = vector<int>(variables.size(), -1);
    }
    set_partial_state_from_action(pre_eff_state, op);
    return pre_eff_state;
}

unique_ptr<Permutation> SymmetricRootTask::get_permutation_for_operator(const ExplicitOperator &op) const {
    vector<int> pre_eff_state(get_operator_post_condition(op));
    unique_ptr<Permutation> perm;
    if (compute_perfect_canonical) {
        perm = make_unique<Permutation>(group->get_perfect_canonical_permutation(pre_eff_state));
    } else {
        perm = make_unique<Permutation>(group->get_canonical_permutation(pre_eff_state));
    }
    return perm;
}

void SymmetricRootTask::set_symmetry_effects_of_operator(const ExplicitOperator &orig_op,
                                                         ExplicitOperator &new_op,
                                                         const unique_ptr<Permutation> &perm) const {
    // TODO implement this
    //  try several variants:
    //  1) get canonical projected onto orig_op's variables => done
    //  2) fill in remaining variables with initial-state values => done
    //  3) fill in remaining variables with random values (same for every action?) => done
    //  4) do a form of context splitting and introduce multiple copies of every action with different symmetries => done
    //  5) incorporate mutex information in context splitting and prune operators with mutex precondition => done
    //  6) do clever context splitting by not doing the full product if variable subsets are not affected by all generators => done
    //  7) decouple computation of canonical permutation if permutations affect different components of variables => done
    //  8) allow for partial splitting to avoid excessive overhead => done

    assert(new_op.effects.empty());

    if (perm->identity()) {
        // no symmetries
        new_op.effects = orig_op.effects;
        return;
    }

    // check for which variables we added effects because of the permutation
    vector<bool> handled_var(original_root_task->get_num_variables(), false);
    vector<int> pre_eff_state(variables.size(), -1);

    if (context_splitting){
        // first we write the new preconditions into the state (defined over all vars)
        // then we write the original effects (effects of new_op are empty at this point)
        set_partial_state_from_action(pre_eff_state, new_op);
        set_partial_state_from_action(pre_eff_state, orig_op);
        assert(skip_unaffected_variables || std::all_of(pre_eff_state.begin(), pre_eff_state.end(), [](int val){return val != -1;}));
    } else {
        set_partial_state_from_action(pre_eff_state, orig_op);
    }

    auto &affected_vars_cycles = perm->affected_vars_cycles;
    for (size_t i = 0; i < affected_vars_cycles.size(); i++) {
        if (affected_vars_cycles[i].size() == 1) {
            int from_var = affected_vars_cycles[i][0];
            int from_val = pre_eff_state[from_var];
            if (from_val != -1) {
                auto [to_var, to_val] = perm->get_new_var_val_by_old_var_val(from_var, from_val);
                // effect is fixed by original precondition or effect
                new_op.effects.push_back({to_var, to_val, {}});
                handled_var[to_var] = true;
                continue;
            }
            assert(skip_unaffected_variables || !context_splitting);
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
            assert(skip_unaffected_variables || !context_splitting);
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
        assert(skip_unaffected_variables || !context_splitting);
        add_conditional_permuted_effects(new_op, *perm, last_var, variables[last_var].domain_size);
        auto [to_var, _] = perm->get_new_var_val_by_old_var_val(last_var, 0);
        handled_var[to_var] = true;
    }

    // copy original effects of variables not affected by the permutation
    for (const auto &eff: orig_op.effects) {
        if (!handled_var[eff.fact.var]) {
            new_op.effects.push_back(eff);
        }
    }
}

void SymmetricRootTask::create_operator(int op_id) {
    const auto &op = original_root_task->operators[op_id];

    ExplicitOperator new_op(op.cost, op.name, op.preconditions); // TODO avoid copying preconditions
    const auto &perm = get_permutation_for_operator(op);
    set_symmetry_effects_of_operator(op, new_op, perm);
    assert(!new_op.effects.empty());

    operators.push_back(new_op);
}

void SymmetricRootTask::add_context_split_cond_effs_recursive(
        size_t var_id,
        std::vector<FactPair> &cond_eff_preconditions,
        const std::vector<int> &component_split_vars,
        const std::vector<ExplicitEffect> &component_effects,
        ExplicitOperator &new_op) {

    if (var_id == component_split_vars.size()){
        ExplicitOperator tmp_op(new_op.cost, new_op.name, cond_eff_preconditions);
        tmp_op.effects = component_effects; // need the effects to compute the permutation
        auto perm = get_permutation_for_operator(tmp_op);

        if (perm->identity() && component_split_vars.empty()){
            std::copy(component_effects.begin(),
                      component_effects.end(),
                      back_inserter(new_op.effects));
            return;
        }

        int num_original_pre = static_cast<int>(cond_eff_preconditions.size() - component_split_vars.size());

        vector<bool> is_affected_var(RootTask::get_num_variables(), false);
        vector<bool> is_eff_var(RootTask::get_num_variables(), false);
        // apply permutation to original effects
        for (const auto &eff : component_effects) {
            is_eff_var[eff.fact.var] = true;
            auto [new_var, new_val] = perm->get_new_var_val_by_old_var_val(eff.fact.var, eff.fact.value);
            assert(!is_affected_var[new_var]);
            is_affected_var[new_var] = true;
            new_op.effects.emplace_back(new_var,
                                        new_val,
                                        vector<FactPair>(cond_eff_preconditions.begin() + num_original_pre,
                                                         cond_eff_preconditions.end()));
        }
        // apply permutation to split variables
        for (size_t i = 0; i < component_split_vars.size(); ++i) {
            int var = component_split_vars[i];
            assert(var == cond_eff_preconditions[num_original_pre + i].var);
            int val = cond_eff_preconditions[num_original_pre + i].value;
            auto [new_var, new_val] = perm->get_new_var_val_by_old_var_val(var, val);
            assert(!is_affected_var[new_var]);
            is_affected_var[new_var] = true;
            if (new_var != var || new_val != val) {
                new_op.effects.emplace_back(new_var,
                                            new_val,
                                            vector<FactPair>(cond_eff_preconditions.begin() + num_original_pre,
                                                             cond_eff_preconditions.end()));
            }
        }
        // apply permutation to prevail conditions
        for (int i = 0; i < num_original_pre; ++i){
            int var = cond_eff_preconditions[i].var;
            if (is_eff_var[var]){
                continue;
            }
            if (num_original_pre < static_cast<int>(cond_eff_preconditions.size())) {
                int val = cond_eff_preconditions[i].value;
                auto [new_var, new_val] = perm->get_new_var_val_by_old_var_val(var, val);
                assert(!is_affected_var[new_var]);
                is_affected_var[new_var] = true;
                if (new_var != var || new_val != val) {
                    new_op.effects.emplace_back(new_var,
                                                new_val,
                                                vector<FactPair>(cond_eff_preconditions.begin() + num_original_pre,
                                                                 cond_eff_preconditions.end()));
                }
            }
        }
        // apply permutation to variables affected by the computation that have not been handled above
        for (int var : perm->vars_affected){
            if (is_affected_var[var]) {
                continue;
            }
            int from_var = -1;
            for (int v_tmp : perm->vars_affected){
                auto [to_var, to_val] = perm->get_new_var_val_by_old_var_val(v_tmp, 0);
                if (to_var == var){
                    from_var = v_tmp;
                    break;
                }
            }
            assert(from_var != -1);
            size_t size_before = new_op.effects.size();
            add_conditional_permuted_effects(new_op, *perm, from_var, variables[from_var].domain_size);
            for (size_t i = size_before; i < new_op.effects.size(); ++i){
                std::copy(cond_eff_preconditions.begin() + num_original_pre,
                          cond_eff_preconditions.end(),
                          back_inserter(new_op.effects[i].conditions));
            }
        }
        return;
    }

    int var = component_split_vars[var_id];
    cond_eff_preconditions.emplace_back(var, 0);
    for (int val = 0; val < RootTask::get_variable_domain_size(var); ++val){
        if (skip_mutex_preconditions) {
            bool is_mutex_pre = false;
            FactPair new_pre(var, val);
            for (size_t i = 0; i < cond_eff_preconditions.size() - 1; ++i) {
                const auto &pre = cond_eff_preconditions[i];
                if (RootTask::are_facts_mutex(pre, new_pre)) {
                    is_mutex_pre = true;
                    break;
                }
            }
            if (is_mutex_pre) {
                continue;
            }
        }

        cond_eff_preconditions.back().value = val;
        add_context_split_cond_effs_recursive(var_id + 1,
                                              cond_eff_preconditions,
                                              component_split_vars,
                                              component_effects,
                                              new_op);
    }
    cond_eff_preconditions.pop_back();
}

void SymmetricRootTask::create_operators_context_split_recursive(size_t var_id,
                                                                 vector<FactPair> &precondition,
                                                                 const vector<int> &outside_post_vars,
                                                                 const int original_op_id) {

    if (var_id == outside_post_vars.size()){
        const auto &original_op = original_root_task->operators[original_op_id];
        ExplicitOperator new_op(original_op.cost, original_op.name, precondition);

        new_op.effects = original_op.effects; // need the effects to compute the permutation
        const auto &perm = get_permutation_for_operator(new_op);
        new_op.effects.clear(); // set_symmetry_effects_of_operator expects empty effects
        set_symmetry_effects_of_operator(original_op, new_op, perm);
        assert(!new_op.effects.empty());

        new_op_id_to_original_op_id.push_back(original_op_id);
        operators.push_back(new_op);
        return;
    }

    int var = outside_post_vars[var_id];
    precondition.emplace_back(var, 0);
    for (int val = 0; val < RootTask::get_variable_domain_size(var); ++val){
        if (skip_mutex_preconditions) {
            bool is_mutex_pre = false;
            FactPair new_pre(var, val);
            for (size_t i = 0; i < precondition.size() - 1; ++i) {
                const auto &pre = precondition[i];
                if (RootTask::are_facts_mutex(pre, new_pre)) {
                    is_mutex_pre = true;
                    break;
                }
            }
            if (is_mutex_pre) {
                continue;
            }
        }

        precondition.back().value = val;
        create_operators_context_split_recursive(var_id + 1,
                                                 precondition,
                                                 outside_post_vars,
                                                 original_op_id);
    }
    precondition.pop_back();
}

void SymmetricRootTask::create_operators_context_split(int op_id) {
    const auto &original_op = original_root_task->operators[op_id];

    vector<int> split_vars(get_split_variables(original_op));

    vector<FactPair> precondition = original_op.preconditions;
    create_operators_context_split_recursive(0,
                                             precondition,
                                             split_vars,
                                             op_id);
}

void SymmetricRootTask::create_operators_context_split_decoupled(int op_id) {
    // TODO: instead of doing this per operator, group them by "action schema" which share this precomputation
    const auto &original_op = original_root_task->operators[op_id];

    vector<int> split_vars(get_split_variables(original_op));

    vector<bool> is_relevant_var(get_num_variables(), false); // var is either split var or in post condition
    for (int var : split_vars){
        is_relevant_var[var] = true;
    }
    for (const auto &pre : original_op.preconditions){
        is_relevant_var[pre.var] = true;
    }
    for (const auto &eff : original_op.effects){
        is_relevant_var[eff.fact.var] = true;
    }

    const vector<vector<int>> &permutation_components = group->get_permutation_components();
    // TODO precompute this?
    vector<vector<int>> affected_components;
    for (const auto &c : permutation_components){
        for (int var : c) {
            if (is_relevant_var[var]){
                affected_components.push_back(c);
                break;
            }
        }
    }

    vector<vector<FactPair>> preconditions_by_component(affected_components.size());
    vector<vector<int>> split_vars_by_component(affected_components.size());
    vector<vector<ExplicitEffect>> effects_by_component(affected_components.size());

    vector<int> var_to_components(get_num_variables(), -1);
    for (size_t c = 0; c < affected_components.size(); ++c){
        for (int var : affected_components[c]){
            var_to_components[var] = static_cast<int>(c);
        }
    }

    ExplicitOperator new_op(original_op.cost, original_op.name, original_op.preconditions);

    for (const auto &pre : original_op.preconditions){
        int c = var_to_components[pre.var];
        if (c != -1) {
            // precondition variable is touched by some permutation
            preconditions_by_component[c].push_back(pre);
        }
    }
    for (const auto &eff : original_op.effects){
        int c = var_to_components[eff.fact.var];
        if (c == -1){
            // effect variable is not touched by any permutation
            new_op.effects.push_back(eff);
        } else {
            effects_by_component[c].push_back(eff);
        }
    }
    for (int var : split_vars){
        int c = var_to_components[var];
        assert(c != -1);
        split_vars_by_component[c].push_back(var);
    }

    for (size_t comp = 0; comp < affected_components.size(); ++comp){
        add_context_split_cond_effs_recursive(0,
                                              preconditions_by_component[comp],
                                              split_vars_by_component[comp],
                                              effects_by_component[comp],
                                              new_op);
    }
    operators.push_back(new_op);
}

void SymmetricRootTask::create_operators() {
    for (int op_id = 0; op_id < static_cast<int>(original_root_task->operators.size()); ++op_id) {
        if (context_splitting) {
            if (decoupled_splitting){
                create_operators_context_split_decoupled(op_id);
            } else {
                create_operators_context_split(op_id);
            }
        } else {
            create_operator(op_id);
        }
    }
}

void SymmetricRootTask::release_memory() {
    // TODO: release memory
}

void SymmetricRootTask::dump() const {
    task_properties::dump_task(TaskProxy(*this), true, true);
}


class SymmetricRootTaskFeature : public plugins::TypedFeature<AbstractTask, SymmetricRootTask> {
public:
    SymmetricRootTaskFeature() : TypedFeature("symmetry") {
        document_title("Symmetric task");
        document_synopsis(
                "A symmetry transformation of the root task.");

        add_option<shared_ptr<Group>>(
                "symmetries",
                "method that computes symmetries",
                "structural_symmetries()");
        add_option<EmptyValueStrategy>(
                "empty_value_strategy",
                "How to treat variables not mentioned in operator precondition and effect",
                "none");
        add_option<bool>(
                "compute_perfect_canonical",
                "Computes the perfect canonical for each orbit.",
                "false");
        add_option<bool>(
                "skip_mutex_preconditions",
                "For empty_value_strategy=split_context, do not generated operators with mutex precondition.",
                "true");
        add_option<bool>(
                "skip_unaffected_variables",
                "For empty_value_strategy=split_context, skip variables not affected by any permutation.",
                "true");
        add_option<bool>(
                "skip_unaffected_variables_relevant_permutations",
                "For empty_value_strategy=split_context, skip variables not affected by any permutation that "
                "touches the postcondition of the respective operator.",
                "true");
        add_option<bool>(
                "context_splitting",
                "Do context splitting over partial states that do not intersection with operator post-condition"
                "to compute symmetrie.",
                "true");
        add_option<bool>(
                "decoupled_splitting",
                "For empty_value_strategy=split_context, if split variables are in different permutation components, "
                "then the enumeration of partial states is done independently for each component.",
                "true");
        add_option<int>(
                "max_number_contexts_per_operator",
                "For empty_value_strategy=split_context, this limits the number of contexts that are generated"
                "per operator, i.e., the number of operator copies, respectively conditional effects (for decoupled_splitting=true).",
                "infinity");
        add_option<bool>(
                "write_sas_file",
                "Writes the decoupled task to dec_output.sas and terminates.",
                "false");
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
        {"random", "use random values for variables not set in precondition or effect"},
        {"goal", "use goal values for variables not set in precondition or effect"},
        {"init_goal", "uses goal values where defined otherwise initial state values"},
        {"split_context", "enumerate all partial states to fill up the post condition"}
    });
}

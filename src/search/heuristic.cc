#include "heuristic.h"

#include "evaluation_context.h"
#include "evaluation_result.h"

#include "plugins/plugin.h"
#include "task_utils/task_properties.h"
#include "tasks/cost_adapted_task.h"
#include "tasks/undecoupled_task.h"
#include "tasks/root_task.h"

#include <cassert>
#include <cstdlib>
#include <limits>

using namespace std;

Heuristic::Heuristic(const plugins::Options &opts)
    : Evaluator(opts, true, true, true),
      heuristic_cache(HEntry(NO_VALUE, true)), //TODO: is true really a good idea here?
      cache_evaluator_values(opts.get<bool>("cache_estimates")),
      task(opts.get<shared_ptr<AbstractTask>>("transform")),
      task_proxy(*task),
      g_root_task_proxy(*tasks::g_root_task),
      undecoupled_task(dynamic_pointer_cast<tasks::UndecoupledTask>(task)),
      state_samples(opts.get<int>("state_samples")),
      min_operator_cost(task_properties::get_min_operator_cost(task_proxy)) {
    if (undecoupled_task && cache_evaluator_values) {
        cache_evaluator_values = false;
        log << "Setting cache_estimates=false for undecoupled task transformation!" << endl;
    }
    sampled_states.reserve(state_samples);
}

Heuristic::~Heuristic() {
}

void Heuristic::set_preferred(const OperatorProxy &op) {
    preferred_operators.insert(op.get_ancestor_operator_id(tasks::g_root_task.get()));
}

State Heuristic::convert_ancestor_state(const State &ancestor_state) const {
    return task_proxy.convert_ancestor_state(ancestor_state);
}

void Heuristic::add_options_to_feature(plugins::Feature &feature) {
    add_evaluator_options_to_feature(feature);
    feature.add_option<shared_ptr<AbstractTask>>(
        "transform",
        "Optional task transformation for the heuristic."
        " Currently, adapt_costs() and no_transform() are available.",
        "no_transform()");
    feature.add_option<bool>("cache_estimates", "cache heuristic estimates", "true");
    feature.add_option<int>("state_samples", "#States sampled for undecoupled heuristic", "1");
}

EvaluationResult Heuristic::compute_result(EvaluationContext &eval_context) {
    EvaluationResult result;

    assert(preferred_operators.empty());

    const State &state = eval_context.get_state();
    bool calculate_preferred = eval_context.get_calculate_preferred();

    int heuristic = NO_VALUE;

    if (!calculate_preferred && cache_evaluator_values &&
        heuristic_cache[state].h != NO_VALUE && !heuristic_cache[state].dirty) {
        heuristic = heuristic_cache[state].h;
        result.set_count_evaluation(false);
    } else {
        if (!undecoupled_task) {
            heuristic = compute_heuristic(state);
            if (cache_evaluator_values) {
                heuristic_cache[state] = HEntry(heuristic, false);
            }
        } else {
            heuristic = numeric_limits<int>::max();
            sampled_states.clear();
            undecoupled_task->get_sampled_states(state, state_samples, sampled_states);
            for (const State &sampled_state : sampled_states) {
                // We encountered a goal state and set the h-value to 0
                if (task_properties::is_goal_state(task_proxy, sampled_state)) {
                    heuristic = 0;
                    break;
                }

                if (our_h_cache.count(sampled_state) == 0) {
                    our_h_cache[sampled_state] = compute_heuristic(sampled_state);
                }
                if (our_h_cache[sampled_state] > DEAD_END) {
                    heuristic = min(our_h_cache[sampled_state], heuristic);
                }
            }

            // Potentially a deadend but not sure
            if (heuristic == numeric_limits<int>::max()) {
                // Blind heuristic
                if (task_properties::is_goal_state(g_root_task_proxy, state)) {
                    heuristic = 0;
                } else {
                    heuristic = min_operator_cost;
                }
            }
        }
        result.set_count_evaluation(true);
    }

    assert(heuristic == DEAD_END || heuristic >= 0);

    if (heuristic == DEAD_END) {
        /*
          It is permissible to mark preferred operators for dead-end
          states (thus allowing a heuristic to mark them on-the-fly
          before knowing the final result), but if it turns out we
          have a dead end, we don't want to actually report any
          preferred operators.
        */
        preferred_operators.clear();
        heuristic = EvaluationResult::INFTY;
    }

    // Remove -1 from preffered operators which can happen for landmark heuristic
    if (undecoupled_task) {
        ordered_set::OrderedSet<OperatorID> cleaned_preferred_operators;
        for (OperatorID op_id : preferred_operators) {
            if (op_id.get_index() != -1) {
                cleaned_preferred_operators.insert(op_id);
            }
        }
        preferred_operators = cleaned_preferred_operators;
    }

    assert(!preferred_operators.contains(OperatorID(-1)));

#ifndef NDEBUG
    TaskProxy global_task_proxy = state.get_task();
    OperatorsProxy global_operators = global_task_proxy.get_operators();
    if (heuristic != EvaluationResult::INFTY) {
        for (OperatorID op_id : preferred_operators) {
            assert(task_properties::is_applicable(global_operators[op_id], state));
        }
    }
#endif

    result.set_evaluator_value(heuristic);
    result.set_preferred_operators(preferred_operators.pop_as_vector());
    assert(preferred_operators.empty());

    return result;
}

bool Heuristic::does_cache_estimates() const {
    return cache_evaluator_values;
}

bool Heuristic::is_estimate_cached(const State &state) const {
    return heuristic_cache[state].h != NO_VALUE;
}

int Heuristic::get_cached_estimate(const State &state) const {
    assert(is_estimate_cached(state));
    return heuristic_cache[state].h;
}

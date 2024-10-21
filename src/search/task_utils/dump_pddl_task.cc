#include "dump_pddl_task.h"

#include "task_properties.h"

#include "../tasks/root_task.h"
#include "../utils/logging.h"

#include <math.h>
#include <regex>

using namespace std;

namespace dump_pddl_task {
const string IND = "    ";

static int get_size_of_binary_encoding(int value) {
    return int(ceil(log2(value)));
}

vector<string> get_binary_encoding(const AbstractTask &task, const FactPair &fact) {
    vector<string> encoding;
    bitset<64> bits(fact.value);
    int domain_size = task.get_variable_domain_size(fact.var);
    for (int i = 0; i < get_size_of_binary_encoding(domain_size); ++i) {
        string cur = "(var" + to_string(fact.var) + "_exp" + to_string(i) + ")";
        if (!bits.test(i)) {
            cur = "(not " + cur + ")";
        }
        encoding.push_back(cur);
    }

    return encoding;
}

void extract_all_preconditions(const AbstractTask &task, int op_no, vector<FactPair> &all_preconditions) {
    for (int pre_ind = 0; pre_ind < task.get_num_operator_preconditions(op_no, false); ++pre_ind) {
        FactPair fact = task.get_operator_precondition(op_no, pre_ind, false);
        all_preconditions.push_back(fact);
    }
}

void extract_all_effects_with_conditions(const AbstractTask &task,
                                         int op_no,
                                         vector<FactPair> &all_effects,
                                         vector<vector<FactPair>> &all_effects_cond) {
    for (int eff_ind = 0; eff_ind < task.get_num_operator_effects(op_no, false); ++eff_ind) {
        vector<FactPair> cond;
        for (int cond_ind = 0; cond_ind < task.get_num_operator_effect_conditions(op_no, eff_ind, false); ++cond_ind) {
            FactPair fact = task.get_operator_effect_condition(op_no, eff_ind, cond_ind, false);
            cond.push_back(fact);
        }
        FactPair eff = task.get_operator_effect(op_no, eff_ind, false);
        all_effects.push_back(eff);
        all_effects_cond.push_back(cond);
    }
}

void dump_domain_header(const AbstractTask & /*task*/, ostream &os) {
    os << "(define (domain dec-domain)" << endl;
}

void dump_domain_footer(const AbstractTask & /*task*/, ostream &os) {
    os << ")" << endl;
}

void dump_domain_requirements(const AbstractTask & /*task*/, ostream &os) {
    os << IND << "(:requirements :adl :derived-predicates :action-costs)" << endl;
}

void dump_domain_predicates(const AbstractTask &task, ostream &os) {
    os << IND << "(:predicates" << endl;
    for (int var = 0; var < task.get_num_variables(); ++var) {
        int domain_size = task.get_variable_domain_size(var);
        int next_power_of_two = int(pow(2, ceil(log2(domain_size))));

        for (const string &var_encoding : get_binary_encoding(task, FactPair(var, next_power_of_two - 1))) {
            assert(!var_encoding.starts_with("(not"));
            os << IND << IND << var_encoding << endl;
        }
    }
    os << IND << ")" << endl;
}

void dump_domain_functions(const AbstractTask & /*task*/, ostream &os) {
    os << IND << "(:functions" << endl;
    os << IND << IND << "(total-cost) - number" << endl;
    os << IND << ")" << endl;
}

void dump_domain_axiom(const AbstractTask &task, ostream &os, int ax_no) {
    assert(task.get_num_operator_effects(ax_no, true) == 1);
    assert(task.get_num_operator_preconditions(ax_no, true) == 1);
    assert(task.get_operator_precondition(ax_no, 0, true).var ==
           task.get_operator_effect(ax_no, 0, true).var);
    assert(task.get_operator_precondition(ax_no, 0, true).value !=
           task.get_operator_effect(ax_no, 0, true).value);

    int var = task.get_operator_effect(ax_no, 0, true).var;
    int val = task.get_operator_effect(ax_no, 0, true).value;

    os << IND << "(:derived " << flush;
    os << get_binary_encoding(task, FactPair(var, val)).at(0) << endl;
    assert(task.get_variable_domain_size(var) == 2);

    os << IND << IND << "(and" << flush;
    for (int cond_ind = 0; cond_ind < task.get_num_operator_effect_conditions(ax_no, 0, true); ++cond_ind) {
        FactPair fact = task.get_operator_effect_condition(ax_no, 0, cond_ind, true);
        for (const string &binary_var : get_binary_encoding(task, fact)) {
            os << " " << binary_var << flush;
        }
    }
    os << ")" << endl;

    os << IND << ")" << endl;
}

void dump_domain_axioms(const AbstractTask &task, ostream &os) {
    for (int ax = 0; ax < task.get_num_axioms(); ++ax) {
        os << endl;
        dump_domain_axiom(task, os, ax);
    }
}

void dump_domain_operator(const AbstractTask &task, ostream &os, int op_no) {
    vector<FactPair> all_preconditions;
    vector<FactPair> all_effects;
    vector<vector<FactPair>> all_effects_cond;

    extract_all_preconditions(task, op_no, all_preconditions);
    extract_all_effects_with_conditions(task, op_no, all_effects, all_effects_cond);

    assert(all_effects.size() == all_effects_cond.size());

    string name = task.get_operator_name(op_no, false);
    // replace all empty spaces with "-----"
    name = regex_replace(name, regex("\\s"), "-----");

    os << IND << "(:action " << name << endl;
    os << IND << IND << ":parameters  ()" << endl;

    os << IND << IND << ":precondition (and" << endl;
    for (const FactPair &fact : all_preconditions) {
        for (const string &binary_var : get_binary_encoding(task, fact)) {
            os << IND << IND << IND << binary_var << endl;
        }
    }
    os << IND << IND << ")" << endl;

    os << IND << IND << ":effect (and" << endl;
    for (size_t eff_id = 0; eff_id < all_effects.size(); ++eff_id) {
        os << IND << IND << IND << flush;
        if (all_effects_cond[eff_id].size() > 0) {
            os << "(when (and" << flush;
            for (const FactPair &fact : all_effects_cond[eff_id]) {
                for (const string &binary_var : get_binary_encoding(task, fact)) {
                    os << " " << binary_var << flush;
                }
            }
            os << ") (and" << flush;
            for (const string &binary_var : get_binary_encoding(task, all_effects[eff_id])) {
                os << " " << binary_var << flush;
            }
            os << "))" << endl;
        } else {
            for (const string &binary_var : get_binary_encoding(task, all_effects[eff_id])) {
                os << binary_var << " " << flush;
            }
            os << endl;
        }
    }
    os << IND << IND << IND << "(increase (total-cost) " << task.get_operator_cost(op_no, false) << ")" << endl;
    os << IND << IND << ")" << endl;

    os << IND << ")" << endl;
}

void dump_domain_operators(const AbstractTask &task, ostream &os) {
    for (int op = 0; op < task.get_num_operators(); ++op) {
        os << endl;
        dump_domain_operator(task, os, op);
    }
}

void dump_problem_header(const AbstractTask & /*task*/, ostream &os) {
    os << "(define (problem dec-problem)" << endl;
    os << IND << "(:domain dec-domain)" << endl;
}

void dump_problem_footer(const AbstractTask & /*task*/, ostream &os) {
    os << IND << "(:metric minimize (total-cost))" << endl;
    os << ")" << endl;
}

void dump_problem_initial_state(const AbstractTask &task, ostream &os) {
    os << IND << "(:init" << endl;
    for (int var = 0; var < task.get_num_variables(); ++var) {
        if (task.get_variable_axiom_layer(var) == -1) {
            int val = task.get_initial_state_values().at(var);
            for (const string &binary_var : get_binary_encoding(task, FactPair(var, val))) {
                os << IND << IND << binary_var << endl;
            }
        }
    }
    os << IND << IND << "(= (total-cost) 0)" << endl;
    os << IND << ")" << endl;
}

void dump_problem_goal(const AbstractTask &task, ostream &os) {
    os << IND << "(:goal (and" << endl;
    for (int goal = 0; goal < task.get_num_goals(); ++goal) {
        for (const string &binary_var : get_binary_encoding(task, task.get_goal_fact(goal))) {
            os << IND << IND << binary_var << endl;
        }
    }
    os << IND << "))" << endl;
}

void dump_domain_as_PDDL(const AbstractTask &task, ostream &os) {
    dump_domain_header(task, os);
    dump_domain_requirements(task, os);
    dump_domain_predicates(task, os);
    dump_domain_functions(task, os);
    dump_domain_axioms(task, os);
    dump_domain_operators(task, os);
    dump_domain_footer(task, os);
}

void dump_problem_as_PDDL(const AbstractTask &task, ostream &os) {
    dump_problem_header(task, os);
    dump_problem_initial_state(task, os);
    dump_problem_goal(task, os);
    dump_problem_footer(task, os);
}
}

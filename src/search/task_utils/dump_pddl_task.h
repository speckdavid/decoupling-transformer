#ifndef TASK_UTILS_DUMP_PDDL_TASK_H
#define TASK_UTILS_DUMP_PDDL_TASK_H

#include "../abstract_task.h"

#include <ostream>

/*
* We encode each FDR variable v with |D| domain values
* using ceil(log2(|D|)) binary variables.
* Those binary variables function as binary counter.
*/

namespace dump_pddl_task {
std::vector<std::string> get_binary_encoding(const AbstractTask &task, const FactPair &fact, bool skip_negated_facts = false);

void extract_all_preconditions(const AbstractTask &task, int op_no, std::vector<FactPair> &all_preconditions);

void extract_all_effects_with_conditions(const AbstractTask &task,
                                         int op_no,
                                         std::vector<FactPair> &all_effects,
                                         std::vector<std::vector<FactPair>> &all_effects_cond);

void dump_domain_header(const AbstractTask &task, std::ostream &os);
void dump_domain_footer(const AbstractTask &task, std::ostream &os);
void dump_domain_requirements(const AbstractTask &task, std::ostream &os);
void dump_domain_predicates(const AbstractTask &task, std::ostream &os);
void dump_domain_functions(const AbstractTask &task, std::ostream &os);
void dump_domain_axiom(const AbstractTask &task, std::ostream &os, int ax_no);
void dump_domain_axioms(const AbstractTask &task, std::ostream &os);
void dump_domain_operator(const AbstractTask &task, std::ostream &os, int op_no);
void dump_domain_operators(const AbstractTask &task, std::ostream &os);

void dump_problem_header(const AbstractTask &task, std::ostream &os);
void dump_problem_footer(const AbstractTask &task, std::ostream &os);
void dump_problem_header(const AbstractTask &task, std::ostream &os);
void dump_problem_initial_state(const AbstractTask &task, std::ostream &os);
void dump_problem_goal(const AbstractTask &task, std::ostream &os);

extern void dump_domain_as_PDDL(const AbstractTask &task, std::ostream &os);
extern void dump_problem_as_PDDL(const AbstractTask &task, std::ostream &os);
}

#endif

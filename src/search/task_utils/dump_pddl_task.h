#ifndef TASK_UTILS_DUMP_PDDL_TASK_H
#define TASK_UTILS_DUMP_PDDL_TASK_H

#include "../abstract_task.h"

#include <ostream>

namespace dump_pddl_task {
std::string get_var_val_name(const AbstractTask &task, const FactPair &fact);

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



extern void dump_domain_as_PDDL(const AbstractTask &task, std::ostream &os);
extern void dump_problem_as_PDDL(const AbstractTask &task, std::ostream &os);
}

#endif

#ifndef TASK_UTILS_DUMP_SAS_TASK_H
#define TASK_UTILS_DUMP_SAS_TASK_H

#include "../abstract_task.h"

#include <ostream>

namespace dump_sas_task {
void extract_all_preconditions(const AbstractTask &task, int op_no, std::vector<FactPair> &all_preconditions);

void extract_all_effects_with_conditions(const AbstractTask &task,
                                         int op_no,
                                         std::vector<FactPair> &all_effects,
                                         std::vector<std::vector<FactPair>> &all_effects_cond);


void dump_version_as_SAS(const AbstractTask &task, std::ostream &os);

void dump_metric_as_SAS(const AbstractTask &task, std::ostream &os);

void dump_variables_as_SAS(const AbstractTask &task, std::ostream &os);

// Rather slow due to the internal representation of mutexes...
void dump_mutexes_as_SAS(const AbstractTask &task, std::ostream &os);

void dump_initial_state_as_SAS(const AbstractTask &task, std::ostream &os);

void dump_goal_as_SAS(const AbstractTask &task, std::ostream &os);

void dump_operator_pre_post_as_SAS(std::ostream &os, int pre, FactPair eff, const std::vector<FactPair> &eff_cond);

void dump_operator_as_SAS(const AbstractTask &task, std::ostream &os, int op_no);

void dump_operators_as_SAS(const AbstractTask &task, std::ostream &os);

void dump_axiom_as_SAS(const AbstractTask &task, std::ostream &os, int op_no);

void dump_axioms_as_SAS(const AbstractTask &task, std::ostream &os);

extern void dump_as_SAS(const AbstractTask &task, std::ostream &os);
}

#endif

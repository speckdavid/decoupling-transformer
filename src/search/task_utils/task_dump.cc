#include "task_dump.h"

#include "task_properties.h"
#include "../tasks/root_task.h"

using namespace std;

namespace task_dump {
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


void dump_version_as_SAS(const AbstractTask & /*task*/, ostream &os) {
    os << "begin_version" << endl;
    os << tasks::PRE_FILE_VERSION << endl;
    os << "end_version" << endl;
}

void dump_metric_as_SAS(const AbstractTask &task, ostream &os) {
    os << "begin_metric" << endl;
    os << !task_properties::is_unit_cost(TaskProxy(task)) << endl;
    os << "end_metric" << endl;
}

void dump_variables_as_SAS(const AbstractTask &task, ostream &os) {
    os << task.get_num_variables() << endl;
    for (int var = 0; var < task.get_num_variables(); ++var) {
        os << "begin_variable" << endl;
        os << task.get_variable_name(var) << endl;
        os << task.get_variable_axiom_layer(var) << endl;
        os << task.get_variable_domain_size(var) << endl;
        for (int val = 0; val < task.get_variable_domain_size(var); ++val) {
            os << task.get_fact_name(FactPair(var, val)) << endl;
        }
        os << "end_variable" << endl;
    }
}

// Rather slow due to the internal representation of mutexes...
void dump_mutexes_as_SAS(const AbstractTask &task, ostream &os) {
    vector<vector<FactPair>> invariant_groups;

    for (int var = 0; var < task.get_num_variables(); ++var) {
        vector<FactPair> invariant_group;

        for (int val = 0; val < task.get_variable_domain_size(var); ++val) {
            FactPair f(var, val);
            invariant_group.push_back(f);

            for (int var2 = 0; var2 < task.get_num_variables(); ++var2) {
                if (var == var2)
                    continue;
                for (int val2 = 0; val2 < task.get_variable_domain_size(var2); ++val2) {
                    FactPair f2(var2, val2);

                    if (task.are_facts_mutex(f, f2))
                        invariant_group.push_back(f2);
                }
            }
        }

        if ((int)invariant_group.size() > task.get_variable_domain_size(var)) {
            // Keeping non-trivial groups
            invariant_groups.push_back(invariant_group);
        }
    }

    os << invariant_groups.size() << endl;
    for (vector<FactPair> invariant_group : invariant_groups) {
        os << "begin_mutex_group" << endl;
        os << invariant_group.size() << endl;
        for (FactPair fact : invariant_group) {
            os << fact.var << " ";
            os << fact.value << endl;
        }
        os << "end_mutex_group" << endl;
    }
}

void dump_initial_state_as_SAS(const AbstractTask &task, ostream &os) {
    os << "begin_state" << endl;
    for (int val : task.get_initial_state_values())
        os << val << endl;
    os << "end_state" << endl;
}

void dump_goal_as_SAS(const AbstractTask &task, ostream &os) {
    os << "begin_goal" << endl;
    os << task.get_num_goals() << endl;
    for (int i = 0; i < task.get_num_goals(); ++i) {
        FactPair fact = task.get_goal_fact(i);
        os << fact.var << " " << fact.value << endl;
    }
    os << "end_goal" << endl;
}

void dump_operator_pre_post_as_SAS(ostream &os, int pre, FactPair eff, const vector<FactPair> &eff_cond) {
    os << eff_cond.size() << " ";
    for (FactPair cond : eff_cond) {
        os << cond.var << " " << cond.value << " " << endl;
    }
    os << eff.var << " " << pre << " " << eff.value << endl;
}

void dump_operator_as_SAS(const AbstractTask &task, ostream &os, int op_no) {
    vector<FactPair> all_preconditions;
    vector<FactPair> all_effects;
    vector<vector<FactPair>> all_effects_cond;

    extract_all_preconditions(task, op_no, all_preconditions);
    extract_all_effects_with_conditions(task, op_no, all_effects, all_effects_cond);

    vector<FactPair> prevail;
    vector<int> eff_pre_val(all_effects.size(), -1);

    for (FactPair c : all_preconditions) {
        // Checking if in any effect
        bool found = false;
        for (size_t i = 0; i < all_effects.size(); ++i) {
            const FactPair &eff = all_effects[i];
            if (eff.var != c.var)
                continue;
            found = true;
            eff_pre_val[i] = c.value;
        }

        if (!found)
            prevail.push_back(c);
    }

    os << "begin_operator" << endl;
    os << task.get_operator_name(op_no, false) << endl;
    os << prevail.size() << endl;
    for (FactPair cond : prevail) {
        os << cond.var << " " << cond.value << endl;
    }
    os << all_effects.size() << endl;
    for (size_t i = 0; i < all_effects.size(); ++i) {
        dump_operator_pre_post_as_SAS(os, eff_pre_val[i], all_effects[i], all_effects_cond[i]);
    }
    os << task.get_operator_cost(op_no, false) << endl;
    os << "end_operator" << endl;
}

void dump_operators_as_SAS(const AbstractTask &task, ostream &os) {
    os << task.get_num_operators() << endl;
    for (int op_no = 0; op_no < task.get_num_operators(); ++op_no) {
        dump_operator_as_SAS(task, os, op_no);
    }
}

void dump_axiom_as_SAS(const AbstractTask &task, ostream &os, int op_no) {
    vector<FactPair> cond;
    for (int cond_ind = 0; cond_ind < task.get_num_operator_effect_conditions(op_no, 0, true); ++cond_ind) {
        FactPair fact = task.get_operator_effect_condition(op_no, 0, cond_ind, true);
        cond.push_back(fact);
    }

    FactPair eff = task.get_operator_effect(op_no, 0, true);
    int eff_pre_val = -1;

    for (int pre_ind = 0; pre_ind < task.get_num_operator_preconditions(op_no, true); ++pre_ind) {
        FactPair fact = task.get_operator_precondition(op_no, pre_ind, true);
        if (eff.var == fact.var) {
            eff_pre_val = fact.value;
            break;
        }
    }

    os << "begin_rule" << endl;
    dump_operator_pre_post_as_SAS(os, eff_pre_val, eff, cond);
    os << "end_rule" << endl;
}

void dump_axioms_as_SAS(const AbstractTask &task, ostream &os) {
    os << task.get_num_axioms() << endl;
    for (int op_no = 0; op_no < task.get_num_axioms(); ++op_no) {
        dump_axiom_as_SAS(task, os, op_no);
    }
}

void dump_as_SAS(const AbstractTask &task, ostream &os) {
    dump_version_as_SAS(task, os);
    dump_metric_as_SAS(task, os);
    dump_variables_as_SAS(task, os);
    dump_mutexes_as_SAS(task, os);
    dump_initial_state_as_SAS(task, os);
    dump_goal_as_SAS(task, os);
    dump_operators_as_SAS(task, os);
    dump_axioms_as_SAS(task, os);
}
}

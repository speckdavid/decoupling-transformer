#include "decoupled_plan_reconstruction_task.h"

#include "../task_utils/task_properties.h"
#include "../plugins/plugin.h"
#include "../decoupling/manual_factoring.h"

#include <iostream>
#include <fstream>
#include <regex>
#include <string>

using namespace std;

namespace tasks {
static string read_file_to_string(const string &filename) {
    ifstream file(filename);
    string content;

    if (file.is_open()) {
        getline(file, content);
        file.close();
    } else {
        cerr << "Unable to open file: " << filename << endl;
    }

    return content;
}


DecoupledPlanReconstructionTask::DecoupledPlanReconstructionTask(const plugins::Options & /*options*/)
    : RootTask(), factoring(nullptr) {
    TaskProxy original_task_proxy(*tasks::g_root_task);
    task_properties::verify_no_axioms(original_task_proxy);
    task_properties::verify_no_conditional_effects(original_task_proxy);

    vector<vector<int>> leaves;
    string leave_str = read_file_to_string("factoring.txt");
    assert(leave_str.size() > 2);
    leave_str = leave_str.substr(1, leave_str.length() - 2);
    leave_str = regex_replace(leave_str, regex("\\],"), "];");
    for (auto leaf_str : utils::split(leave_str, ";")) {
        leaves.push_back(vector<int>());
        leaf_str = leaf_str.substr(1, leaf_str.length() - 2);
        for (auto var_str : utils::split(leaf_str, ",")) {
            leaves.back().push_back(stoi(var_str));
        }
    }

    plugins::Options opts;
    opts.set("verbosity", utils::Verbosity::NORMAL);
    opts.set("min_number_leaves", 1);
    opts.set("max_leaf_size", numeric_limits<int>::max());
    opts.set("factoring_time_limit", numeric_limits<int>::max());
    opts.set("prune_fork_leaf_state_spaces", false);
    opts.set("leaves", leaves);

    factoring = make_shared<decoupling::ManualFactoring>(opts);

    factoring->compute_factoring();

    cout << "TODO!!!" << endl;
    exit(0);
}

class DecoupledPlanReconstructionTaskFeature : public plugins::TypedFeature<AbstractTask, DecoupledPlanReconstructionTask> {
public:
    DecoupledPlanReconstructionTaskFeature() : TypedFeature("decoupled_plan_reconstruction") {
        document_title("Decoupled plan reconstruction task");
        document_synopsis("Task transformation that reconstructs a plan (sas_plan) for a given problem (output.sas) and factoring (factoring.txt)");
    }

    virtual shared_ptr<DecoupledPlanReconstructionTask> create_component(const plugins::Options &options, const utils::Context &) const override {
        return make_shared<DecoupledPlanReconstructionTask>(options);
    }
};

static plugins::FeaturePlugin<DecoupledPlanReconstructionTaskFeature> _plugin;
}

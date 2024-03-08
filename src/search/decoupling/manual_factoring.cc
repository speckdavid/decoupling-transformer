#include "manual_factoring.h"

#include "../plugins/plugin.h"
#include "../task_utils/causal_graph.h"


using namespace std;


namespace decoupling {
ManualFactoring::ManualFactoring(const plugins::Options &opts) :
    Factoring(opts) {
    min_number_leaves = 1;
    leaves = opts.get_list<vector<int>>("leaves");
    if (log.is_at_least_normal()) {
        log << "Using Manual factoring method: " << leaves << endl;
    }
}

void ManualFactoring::compute_factoring_() {
    compute_action_schemas();

    if (action_schemas.empty()) {
        // mostly for trivially unsolvable task from translator?
        log << "ERROR: No action schemas." << endl;
        return;
    }

    if (!is_valid_factoring(leaves)) {
        utils::g_log << "Leaves: " << leaves << endl;
        cerr << "Not a valid factoring!" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }
}

bool ManualFactoring::is_valid_factoring(const vector<vector<int>> &leaves) const {
    unordered_set<int> seen_vars;
    int summed_size = 0;

    for (const vector<int> &leaf : leaves) {
        summed_size += leaf.size();
        for (int var : leaf) {
            if (var < 0 || var >= (int)task_proxy.get_variables().size()) {
                return false;
            }
            seen_vars.insert(var);
        }
    }

    return (int)seen_vars.size() == summed_size;
}


class ManualFactoringFeature : public plugins::TypedFeature<Factoring, ManualFactoring> {
public:
    ManualFactoringFeature() : TypedFeature("manual") {
        document_title("Manual factoring");

        Factoring::add_options_to_feature(*this);
        this->add_list_option<vector<int>>("leaves", "Leaves of the factoring.");
    }
};

static plugins::FeaturePlugin<ManualFactoringFeature> _plugin;
}

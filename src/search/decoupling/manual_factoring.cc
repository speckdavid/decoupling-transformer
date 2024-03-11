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
    // nothing to do
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

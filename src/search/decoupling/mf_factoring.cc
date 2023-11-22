#include "mf_factoring.h"

#include "../plugins/plugin.h"
#include "../task_utils/causal_graph.h"

#include <algorithm>
#include <iostream>
#include <queue>


using namespace std;


namespace decoupling {

MFFactoring::MFFactoring(const plugins::Options &opts) : Factoring(opts) {
    min_number_leaves = 1;
    if (log.is_at_least_normal()) {
        log << "Using Miura & Fukunaga factoring method." << endl;
    }
}

inline bool subseteq(const vector<int> &a, const vector<int> &b) {
    if (a.size() > b.size()){
        return false;
    }
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()){
        if (a[i] < b[j]){
            return false;
        } else if (a[i] > b[j]){
            ++j;
        } else if (a[i] == b[j]){
            ++i;
            ++j;
        }
    }
    if (i < a.size()){
        return false;
    }
    return true;
}

inline vector<int> intersection(const vector<int> &a, const vector<int> &b) {
    vector<int> inters;
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()){
        if (a[i] < b[j]){
            ++i;
        } else if (a[i] > b[j]){
            ++j;
        } else if (a[i] == b[j]){
            inters.push_back(a[i]);
            ++i;
            ++j;
        }
    }
    return inters;
}

void MFFactoring::compute_factoring_() {

    compute_action_schemas();

    if (action_schemas.empty()){
        // mostly for trivially unsolvable task from translator?
        log << "ERROR: No action schemas." << endl;
        return;
    }

    priority_queue<vector<int>> candidates;

    // initial candidates are the precondition variables of all actions
    for (const auto &as : action_schemas){
        candidates.push(as.pre_vars);
    }

    utils::HashSet<vector<int>> closed_candidates;

    vector<int> best_candidate;
    while (!candidates.empty()) {
        auto &c = candidates.top();
        if (closed_candidates.count(c) > 0){
            candidates.pop();
            continue;
        }
        closed_candidates.insert(c);
        bool is_valid = true;
        for (const auto &as: action_schemas) {
            if (!subseteq(c, as.pre_vars) && !subseteq(as.eff_vars, c)) {
                auto pre_intersect = intersection(as.pre_vars, c);
                auto eff_intersect = intersection(as.eff_vars, c);

                if (!pre_intersect.empty() && closed_candidates.count(c) == 0) {
                    candidates.push(pre_intersect);
                }
                if (!eff_intersect.empty() && closed_candidates.count(c) == 0) {
                    candidates.push(eff_intersect);
                }

                is_valid = false;
                break;
            }
        }
        if (is_valid && c.size() > best_candidate.size()){
            best_candidate = c;
        }
        candidates.pop();
    }

    if (best_candidate.empty()){
        log << "No valid decomposition possible." << endl;
        return;
    }
    leaves = vector<vector<int>>(1, best_candidate);
}

class MFFactoringFeature : public plugins::TypedFeature<Factoring, MFFactoring> {
public:
    MFFactoringFeature() : TypedFeature("mf") {
        document_title("Miura & Fukunaga factoring");

        Factoring::add_options_to_feature(*this);
    }
};

static plugins::FeaturePlugin<MFFactoringFeature> _plugin;
}

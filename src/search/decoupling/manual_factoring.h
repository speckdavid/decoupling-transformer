#ifndef DECOUPLING_MANUAL_FACTORING
#define DECOUPLING_MANUAL_FACTORING

#include "factoring.h"

#include "../plugins/options.h"

namespace decoupling {
// manual method to specify factorings by hand
class ManualFactoring : public decoupling::Factoring {
    void compute_factoring_() override;

public:

    explicit ManualFactoring(const plugins::Options &opts);
};
}

#endif

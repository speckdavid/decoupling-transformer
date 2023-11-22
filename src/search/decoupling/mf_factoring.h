#ifndef DECOUPLING_MF_FACTORING
#define DECOUPLING_MF_FACTORING

#include "factoring.h"

#include "../plugins/options.h"

namespace decoupling {

// "factoring" method from Miura & Fukunaga (ICAPS 2017)
class MFFactoring : public decoupling::Factoring {

    void compute_factoring_() override;

public:

    explicit MFFactoring(const plugins::Options &opts);

};
}

#endif

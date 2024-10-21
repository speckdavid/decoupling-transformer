#include "leaf_state_id.h"

#include <limits>
#include <ostream>

using namespace std;

namespace decoupling {
const LeafStateHash LeafStateHash::MAX(numeric_limits<size_t>::max());

const FactorID FactorID::CENTER(numeric_limits<unsigned short>::max());

const LeafStateID LeafStateID::no_state(LeafStateHash::MAX, FactorID::CENTER);

ostream &operator<<(ostream &os, const LeafStateID &id) {
    os << "#" << id.hash() << " factor " << id.get_factor();
    return os;
}
}

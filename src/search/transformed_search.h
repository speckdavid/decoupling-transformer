#ifndef TRANSFORMED_SEARCH_H
#define TRANSFORMED_SEARCH_H

#include "search_algorithm.h"

#include <memory>
#include <vector>

namespace transformed_search {
class TransformedSearch : public SearchAlgorithm {
public:
    TransformedSearch(const plugins::Options &opts);
    virtual ~TransformedSearch();
};
}

#endif

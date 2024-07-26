#ifndef DECOUPLING_LEAF_STATE_ID_H
#define DECOUPLING_LEAF_STATE_ID_H

#include <iostream>

// For documentation on classes relevant to storing and working with registered
// states see the file state_registry.h.

namespace decoupling {
struct LeafStateHash {
private:
    size_t val;
public:
    explicit LeafStateHash(size_t value) : val(value) {}
    operator size_t() const {
        return val;
    }
    void operator++() {
        ++val;
    }
    bool operator==(const LeafStateHash other) const {
        return val == other.val;
    }
    static const LeafStateHash MAX;
};

//template <>
//struct std::hash<LeafStateHash> {
//    std::size_t operator()(const LeafStateHash & key) const {
//        return key;
//    }
//};

struct FactorID {
private:
    unsigned short val;
public:
    explicit FactorID(unsigned short value) : val(value) {}
    operator unsigned short() const {
        return val;
    }
    void operator++() {
        ++val;
    }
    bool operator==(const FactorID other) const {
        return val == other.val;
    }
    static const FactorID CENTER;
};

//template <>
//struct std::hash<FactorID> {
//    std::size_t operator()(const FactorID & key) const {
//        return key;
//    }
//};

class LeafStateID {
    friend class StateRegistry;
    friend std::ostream &operator<<(std::ostream &os, LeafStateID id);

    LeafStateHash value;

    FactorID factor;

    // No implementation to prevent default construction
    LeafStateID();

public:

    explicit LeafStateID(LeafStateHash id, FactorID factor)
        : value(id), factor(factor) {
    }

    bool operator==(const LeafStateID &other) const {
        return value == other.value;
    }

    bool operator!=(const LeafStateID &other) const {
        return !(*this == other);
    }

    bool operator<(const LeafStateID &other) const {
        return value < other.value;
    }

    LeafStateHash hash() const {
        return value;
    }

    FactorID get_factor() const {
        return factor;
    }


    static const LeafStateID no_state;
};

std::ostream &operator<<(std::ostream &os, const LeafStateID &id);
}

#endif

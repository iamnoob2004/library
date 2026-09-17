#pragma once

template<typename T>
struct monoid_min{
    using value_type=T;
    using S=value_type;

    static constexpr S op(const S &a, const S &b) {
        if(a<b) return a;
        return b;
    }

    static constexpr S id(){return inf<T>;}
};
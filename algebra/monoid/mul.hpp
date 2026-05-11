#pragma once

template<typename T>
struct monoid_mul{
    using value_type=T;

    static constexpr T op(const T &a, const T &b) {return a*b;}
    static constexpr T id(){return T(1);}

    static constexpr T inv(const T &a) {return T(1)/a;}
};
#pragma once

template<typename T>
struct monoid_add{
    using value_type=T;

    static constexpr T op(const T &a, const T &b) {return a+b;}
    static constexpr T id(){return T(0);}

    static constexpr T inv(const T &a) {return -a;}
};
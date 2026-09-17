#pragma once

#include "algebra/monoid/add.hpp"
#include "algebra/monoid/max.hpp"

template<typename T>
struct monoid_action_max_add{
    using monoid_S=monoid_max<T>;
    using monoid_F=monoid_add<T>;

    using S=typename monoid_S::value_type;
    using F=typename monoid_F::value_type;

    static constexpr S act(const F &f, const S& s){
        return s+f;
    }
};
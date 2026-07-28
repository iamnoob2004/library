#pragma once

#include "algebra/monoid/add.hpp"
#include "algebra/monoid/mincnt.hpp"

template<typename T1, typename T2>
struct monoid_action_mincnt_add{
    using monoid_S=monoid_mincnt<T1,T2>;
    using monoid_F=monoid_add<T1>;

    using S=typename monoid_S::value_type;
    using F=typename monoid_F::value_type;

    static constexpr S act(const F &f, const S& s){
        return {s.first+f,s.second};
    }
};
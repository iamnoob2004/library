#pragma once

template<typename T1, typename T2>
struct monoid_mincnt{
    using value_type=pair<T1,T2>;
    using S=value_type;

    static constexpr S op(const S &a, const S &b) {
        if(a.first<b.first) return a;
        if(b.first<a.first) return b;
        return {a.first,a.second+b.second};
    }

    static constexpr S id(){return {inf<T1>,T2(0)};}
};
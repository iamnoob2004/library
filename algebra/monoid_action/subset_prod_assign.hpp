#pragma once

#include "algebra/monoid/add_n.hpp"
#include "algebra/monoid/assign_n.hpp"

template<typename T, int N, int none_val>
struct monoid_action_subset_prod_assign{
    using value_type=T;
    using monoid_S=monoid_add_n<T,(1<<N)>;
    using monoid_F=monoid_assign_n<T,N,none_val>;

    using S=typename monoid_S::S;
    using F=typename monoid_F::S;

    static S act(const F &f, const S& s){
        S res=s;
        int msk=0;
        for(int i=0; i<N; ++i){
            if(f[i]!=T(none_val)){
                msk|=1<<i;
            }
        }
        for(int i=1; i<(1<<N); ++i){
            int msk2=i&msk;
            if(msk2){
                int j=__builtin_ctz(msk2);
                res[i]=res[i-(1<<j)]*f[j];
            }
        }
        return res;
    }
};
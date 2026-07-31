#pragma once

template<typename T, int N>
struct monoid_add_n{
    using value_type=T;
    using S=array<T,N>;

    static S op(const S &a, const S &b) {
        S res=a;
        for(int i=0; i<N; ++i){
            res[i]=a[i]+b[i];
        }
        return res;
    }

    static S id(){S res; for(int i=0; i<N; ++i) res[i]=T(0); return res;}
};
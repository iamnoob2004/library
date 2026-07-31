#pragma once

template<typename T, int N, int none_val>
struct monoid_assign_n{
    using value_type=T;
    using S=array<T,N>;

    static S op(const S &a, const S &b) {
        S res=a;
        for(int i=0; i<N; ++i){
            if(a[i]==T(none_val)) res[i]=b[i];
        }
        return res;
    }

    static S id(){S res; for(int i=0; i<N; ++i) res[i]=T(none_val); return res;}
};
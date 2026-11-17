#pragma once

#include "nt/basic.hpp"

inline int primitive_root_prime(int p){
    if(p==2) return 1;
    int m=p-1;
    vector<int> f;
    for(int i=2; i*i<=m; ++i){
        if(m%i==0){
            f.push_back(i);
            while(m%i==0) m/=i;
        }
    }
    if(m) f.push_back(m);
    for(int g=2; ; ++g){
        bool fail=false;
        for(int q: f){
            if(modpow(g,(p-1)/q,p)==1){
                fail=1;
                break;
            }
        }
        if(!fail) return g;
    }
}
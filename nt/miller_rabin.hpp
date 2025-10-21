#pragma once

#include "library/mod/montgomery_modint.hpp"

// https://judge.yosupo.jp/problem/primality_test
bool is_prime(ll x){
    if(x<2) return false;
    static const vector<int> small={2,3,5,7,11,13,17,19};
    for(int p: small){
        if(x==p) return true;
        if(x%p==0) return false;
    }
    if(x<400) return true;
    assert(x<(1ll<<62));
    ll d=x-1;
    int s=0;
    while(d%2==0) d>>=1,s++;
    using mint=montgomery_modint_64<777771449,false>;
    mint::set_mod(x);
    const mint zero(0),one(1),minus_one(x-1);
    auto check=[&](ll _a) -> bool{
        mint a=mint(_a).pow(d);
        if(a==one||a==zero) return true;
        for(int i=0; i<s; ++i,a*=a){
            if(a==one) return false;
            if(a==minus_one) return true;
        }
        return false;
    };
    if(x<(1ll<<32)){
        for(ll a: {2,7,61}){
            if(!check(a)) return false;
        }
    }
    else{
        for(ll a: {2,325,9375,28178,450775,9780504,1795265022}){
            if(!check(a)) return false;
        }
    }
    return true;
}
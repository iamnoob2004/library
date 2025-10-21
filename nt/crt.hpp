#pragma once

#include "library/nt/binary_gcd.hpp"
#include "library/nt/extgcd.hpp"

pair<ll,ll> CRT(ll x1, ll m1, ll x2, ll m2){
    auto [g,val1,val2]=extgcd<ll>(m1,m2);
    if((x2-x1)%g) return {-1,-1};
    m1/=g,m2/=g;
    ll lcm=m1*m2*g;
    __int128 res=((__int128)val1)*(x2-x1)*m1+x1;
    return {(res%lcm+lcm)%lcm,lcm};
}

pair<ll,ll> CRT(vector<ll> vals, vector<ll> mods){
    assert(vals.size()==mods.size());
    ll x=vals[0],m=mods[0];
    for(int i=1; i<(int)vals.size(); ++i){
        auto [nx,nm]=CRT(x,m,vals[i],mods[i]);
        if(nx==-1) return {-1,-1};
        x=nx,m=nm;
    }
    return {x,m};
}
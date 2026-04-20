#pragma once

#include "library/nt/extgcd.hpp"

template<typename mint>
mint garner3(vector<ll> x, vector<ll> m){
    ll r01=modinv<ll>(m[0],m[1]);
    ll r02=modinv<ll>(m[0],m[2]);
    ll r12=modinv<ll>(m[1],m[2]);
    ll res0=x[0];
    ll res1=(x[1]-res0%m[1]+m[1])*r01%m[1];
    ll res2=(x[2]-res0%m[2]+m[2])*r02%m[2];
    res2=(res2-res1%m[2]+m[2])*r12%m[2];
    mint res;
    res+=mint(res0);
    res+=mint(res1)*mint(m[0]);
    res+=mint(res2)*mint(m[0])*mint(m[1]);
    return res;
}
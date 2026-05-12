#pragma once

#include "library/poly/poly_mod_pow.hpp"
#include "library/poly/poly_gcd.hpp"
#include "library/random/basic.hpp"

// find distinct roots
template<typename mint>
vector<mint> poly_find_root(fps<mint> a){
    a.topos();
    assert(!a.empty());
    const int p=mint::get_mod();
    if(p<49){
        vector<mint> res;
        for(int i=0; i<p; ++i){
            if(a.eval(mint(i))==mint(0)){
                res.pb(mint(i));
            }
        }
        return res;
    }
    fps<mint> b={0,1};
    b=poly_mod_pow<mint>(b,p,a);
    b-=fps<mint>({0,1});
    b=poly_gcd<mint>(b,a);

    my_random my_device;

    vector<mint> res;
    auto dfs=[&](auto &dfs, fps<mint> f){
        if((int)f.size()==1) return;
        if((int)f.size()==2){
            res.pb(-f[0]/f[1]);
            return;
        }
        fps<mint> c={my_device.rng(p),1};
        c=poly_mod_pow<mint>(c,(p-1)/2,f);
        c-=fps<mint>({1});
        if((int)c.size()<=1){
            return dfs(dfs,f);
        }
        fps<mint> f1=poly_gcd<mint>(f,c);
        fps<mint> f2=poly_div<mint>(f,f1).first;
        dfs(dfs,f1),dfs(dfs,f2);
    };

    dfs(dfs,b);
    sort(res.begin(),res.end());
    return res;
}
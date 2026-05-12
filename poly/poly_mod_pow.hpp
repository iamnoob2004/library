#pragma once

#include "library/poly/poly_div.hpp"

template<typename mint>
vector<mint> poly_mod_pow(vector<mint> a, ll n, vector<mint> b){
    if(n==0) return {1};
    if(n==1) return a;
    vector<mint> c=poly_mod_pow(a,n/2,b);
    c=convolution<mint>(c,c);
    c=poly_div<mint>(c,b).second;
    if(n&1){
        c=convolution<mint>(c,a);
        c=poly_div<mint>(c,b).second;
    }
    return c;
}
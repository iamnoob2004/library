#pragma once

#include "library/mod/modint_basic.hpp"
#include "library/mod/pow_table.hpp"
#include "library/poly/middle_product.hpp"

// f(x) -> f(x+c)
template<typename mint>
vector<mint> poly_taylor_shift(vector<mint> a, mint c){
    if(c==mint(0)) return a;
    int n=a.size();
    for(int i=0; i<n; ++i) a[i]*=fac<mint>(i);
    vector<mint> b=pow_table(c,n);
    for(int i=0; i<n; ++i) b[i]*=ifac<mint>(i);
    a.resize(n*2-1);
    vector<mint> res=middle_product<mint>(a,b);
    for(int i=0; i<n; ++i) res[i]*=ifac<mint>(i);
    return res;
}
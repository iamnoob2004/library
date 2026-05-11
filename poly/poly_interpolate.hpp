#pragma once

#include "library/algebra/monoid/mul.hpp"
#include "library/ds/sliding_window_aggregation.hpp"
#include "library/mod/modint_basic.hpp"
#include "library/poly/middle_product.hpp"

// given f(0),f(1),...,f(n-1) with deg <= n-1, return f(c)
template<typename mint>
mint lagrange_interpolate_iota(vector<mint> f, mint c){
    int n=f.size();
    if(c.x<n) return f[c.x];
    auto a=f;
    for(int i=0; i<n; ++i){
        a[i]*=ifac<mint>(i)*ifac<mint>(n-1-i);
        if((n-1-i)&1) a[i]=-a[i];
    }
    vector<mint> pref(n+1,1),suff(n+1,1);
    for(int i=0; i<n; ++i) pref[i+1]=pref[i]*(c-mint(i));
    for(int i=n-1; i>=0; --i) suff[i]=suff[i+1]*(c-mint(i));
    mint res;
    for(int i=0; i<n; ++i){
        res+=a[i]*pref[i]*suff[i+1];
    }
    return res;
}

// given f(0),f(1),...,f(n-1) with deg <= n-1, return f(c),f(c+1),...,f(c+m-1)
template<typename mint>
vector<mint> lagrange_interpolate_iota(vector<mint> f, mint c, int m){
    if(m<=49){
        vector<mint> res(m);
        for(int i=0; i<m; ++i) res[i]=lagrange_interpolate_iota<mint>(f,c+mint(i));
        return res;
    }
    int n=f.size();
    auto a=f;
    for(int i=0; i<n; ++i){
        a[i]*=ifac<mint>(i)*ifac<mint>(n-1-i);
        if((n-1-i)&1) a[i]=-a[i];
    }
    vector<mint> b(n+m-1);
    for(int i=0; i<n+m-1; ++i){
        mint x=c+mint(i-(n-1));
        if(x!=mint(0)) b[i]=x.inv();
    }
    reverse(a.begin(),a.end());
    a=middle_product<mint>(b,a);

    sliding_window_aggregation<monoid_mul<mint>> swag;
    for(int i=0; i<n; ++i){
        swag.push(c+mint(i-n));
    }
    vector<mint> res(m);
    for(int i=0; i<m; ++i){
        swag.push(c+mint(i));
        swag.pop();
        mint tot=swag.prod();
        if(tot==mint(0)){
            int j=(c+mint(i)).x;
            assert(j>=0&&j<n);
            res[i]=f[j];
        }
        else{
            res[i]=a[i]*tot;
        }
    }
    return res;
}
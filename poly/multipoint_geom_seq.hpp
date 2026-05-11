#pragma once

#include "library/mod/all_inverse.hpp"
#include "library/poly/middle_product.hpp"

// calculate f(ar^k) for 0 <= k < m
template<typename mint>
vector<mint> multipoint_eval_on_geom_seq(vector<mint> f, mint a, mint r, int m){
    int n=f.size();
    if(m==0) return {};
    auto eval=[&](mint x) -> mint{
        mint res,pw=1;
        for(int i=0; i<n; ++i){
            res+=f[i]*pw;
            pw*=x;
        }
        return res;
    };
    if(r==mint(0)){
        vector<mint> res(m);
        for(int i=1; i<m; ++i) res[i]=f[0];
        res[0]=eval(a);
        return res;
    }
    if(n<49||m<49){
        vector<mint> res(m);
        for(int i=0; i<m; ++i,a*=r) res[i]=eval(a);
        return res;
    }
    {
        mint pw=1;
        for(int i=0; i<n; ++i,pw*=a) f[i]*=pw;
    }
    auto calc=[&](mint q, int k){
        // calc q^binom(i,2) for 0 <= i < k
        vector<mint> res(k);
        mint pw=1;
        res[0]=1;
        for(int i=0; i<k-1; ++i,pw*=q){
            res[i+1]=res[i]*pw;
        }
        return res;
    };
    vector<mint> vec1=calc(r.inv(),max(n,m)),vec2=calc(r,n+m);
    vector<mint> c=f;
    for(int i=0; i<n; ++i){
        c[i]=f[i]*vec1[i];
    }
    c=middle_product<mint>(vec2,c);
    vector<mint> res(m);
    for(int i=0; i<m; ++i){
        res[i]=c[i]*vec1[i];
    }
    return res;
}

// y_i = f(ar^i)
template<typename mint>
vector<mint> multipoint_interpolation_on_geom_seq(vector<mint> y, mint a, mint r){
    int n=y.size();
    if(n==0) return {};
    if(n==1) return {y[0]};
    assert(r!=mint(0));
    mint rinv=r.inv();

    vector<mint> pw(n*2-1),pw2(n*2-1),ipw(n*2-1),ipw2(n*2-1);
    pw[0]=pw2[0]=ipw[0]=ipw2[0]=1;
    for(int i=0; i<n*2-2; ++i){
        pw[i+1]=pw[i]*r;
        pw2[i+1]=pw2[i]*pw[i];
        ipw[i+1]=ipw[i]*rinv;
        ipw2[i+1]=ipw2[i]*ipw[i];
    }

    vector<mint> s(n);
    s[0]=1;
    rep(i,1,n){
        s[i]=s[i-1]*(mint(1)-pw[i]);
    }
    auto sinv=all_inverse<mint>(s);
    s.pb(s[n-1]*(mint(1)-pw[n]));
    for(int i=0; i<n; ++i){
        y[i]*=pw2[n-1-i]*ipw2[n-1]*sinv[i]*sinv[n-1-i]*(i&1?mint(-1):mint(1));
    }

    for(int i=0; i<n; ++i) y[i]*=ipw2[i];
    vector<mint> f=middle_product<mint>(pw2,y);
    for(int i=0; i<n; ++i) f[i]*=ipw2[i];

    vector<mint> g(n);
    g[0]=1;
    for(int i=1; i<n; ++i){
        g[i]=pw2[i]*s[n]*sinv[i]*sinv[n-i]*(i&1?mint(-1):mint(1));
    }
    f=convolution<mint>(f,g);
    f.resize(n);

    reverse(all(f));
    mint ainv=a.inv(),tmp=1;
    for(int i=0; i<n; ++i){
        f[i]*=tmp;
        tmp*=ainv;
    }

    return f;
}

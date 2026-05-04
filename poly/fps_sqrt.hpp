#pragma once

#include "library/poly/convolution.hpp"
#include "library/mod/modint_basic.hpp"
#include "library/poly/fps_basic.hpp"
#include "library/poly/fps_pow.hpp"
#include "library/nt/sqrt_mod.hpp"

// const term = 1: fps_sqrt_1

template<typename mint>
vector<mint> fps_sqrt_1_dense(vector<mint> a){
    assert(!a.empty()&&a[0]==mint(1));
    int n=a.size();
    vector<mint> res={1};
    mint inv2=inv<mint>(2);
    for(int m=1; m<n; m<<=1){
        res.resize(m*2);
        vector<mint> tmp(a.begin(),a.begin()+min(m*2,n));
        tmp=convolution<mint>(tmp,fps_inv<mint>(res));
        tmp.resize(m*2);
        for(int i=0; i<m*2; ++i) res[i]+=tmp[i];
        for(int i=0; i<m*2; ++i) res[i]*=inv2;
    }
    res.resize(n);
    return res;
}

template<typename mint>
vector<mint> fps_sqrt_1_sparse(vector<mint> a){
    return fps_pow_1_sparse<mint>(a,inv<mint>(2));
}

template<typename mint>
vector<mint> fps_sqrt_1(vector<mint> a){
    assert(!a.empty()&&a[0]==mint(1));
    int cnt=nonzero_count(a);
    int thres=mint::can_ntt()?123:1234;
    if(cnt<=thres) return fps_sqrt_1_sparse<mint>(a);
    return fps_sqrt_1_dense<mint>(a);
}

template<typename mint>
vector<mint> fps_sqrt(vector<mint> a){
    if(a.empty()) return {};
    int n=a.size();
    int m=n;
    for(int i=n-1; i>=0; --i){
        if(a[i]!=mint(0)) m=i;
    }
    if(m==n) return a;
    if(m&1) return {};
    if(!is_quadratic_residue<mint>(a[m])) return {};
    int offset=m/2;
    mint c=a[m],cinv=c.inv();
    vector<mint> b(n-offset);
    for(int i=0; i<n-m; ++i) b[i]=a[m+i]*cinv;
    b=fps_sqrt_1<mint>(b);
    vector<mint> res(n);
    c=sqrt_mod<mint>(c);
    for(int i=0; i<n-offset; ++i) res[offset+i]=b[i]*c;
    return res;
}
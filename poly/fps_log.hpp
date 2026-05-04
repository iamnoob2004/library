#pragma once

#include "library/poly/convolution.hpp"
#include "library/mod/modint_basic.hpp"
#include "library/poly/fps_basic.hpp"
#include "library/poly/fps_inv.hpp"

template<typename mint>
vector<mint> fps_log_dense(vector<mint> a){
    int n=(int)a.size();
    assert(n>0&&a[0]==mint(1));
    if(n==1) return {0};
    vector<mint> d=derivative<mint>(a);
    vector<mint> b=a;
    b.pop_back();
    vector<mint> res=convolution<mint>(d,fps_inv<mint>(b));
    res.resize(n-1);
    return integral<mint>(res);
}

template<typename mint>
vector<mint> fps_log_sparse(vector<mint> a){
    int n=(int)a.size();
    assert(n>0&&a[0]==mint(1));
    if(n==1) return {0};
    vector<pair<int,mint>> b; // non-zero
    for(int i=1; i<n; ++i){
        if(a[i]!=mint(0)) b.push_back({i,a[i]});
    }
    vector<mint> c(n-1),res(n);
    for(int i=0; i<n-1; ++i){
        mint val=a[i+1]*mint(i+1);
        for(auto [j,x]: b){
            if(j>i) break;
            val-=c[i-j]*x;
        }
        c[i]=val,res[i+1]=val*inv<mint>(i+1);
    }
    return res;
}

template<typename mint>
vector<mint> fps_log(vector<mint> a){
    int n=(int)a.size();
    assert(n>0&&a[0]!=0);
    int cnt=nonzero_count(a);
    int thres=mint::can_ntt()?234:1234;
    if(cnt<=thres) return fps_log_sparse<mint>(a);
    return fps_log_dense<mint>(a);
}
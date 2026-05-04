#pragma once

#include "library/poly/convolution.hpp"
#include "library/poly/fps_basic.hpp"

template<typename mint>
vector<mint> fps_inv_dense(vector<mint> a){
    int n=(int)a.size();
    assert(n>0&&a[0]!=0);
    vector<mint> res(1,a[0].inv());
    auto b=a;
    for(int m=1; m<n; m<<=1){
        if(n<m*2) b.resize(m*2);
        vector<mint> v1=b,v2=res;
        v1.resize(m*2);
        v1=convolution<mint>(v1,v2),v1.resize(m*2);
        v1=convolution<mint>(v1,v2);
        res.resize(m*2);
        for(int i=0; i<m; ++i) res[i]+=res[i];
        for(int i=0; i<m*2; ++i) res[i]-=v1[i];
    }
    res.resize(n);
    return res;
}

template<typename mint>
vector<mint> fps_inv_sparse(vector<mint> a){
    int n=(int)a.size();
    assert(n>0&&a[0]!=0);
    vector<pair<int,mint>> b; // non-zero
    for(int i=1; i<n; ++i){
        if(a[i]!=mint(0)) b.push_back({i,a[i]});
    }
    mint ainv=a[0].inv();
    vector<mint> res(n);
    res[0]=ainv;
    for(int i=1; i<n; ++i){
        for(auto [j,x]: b){
            if(j>i) break;
            res[i]-=res[i-j]*x;
        }
        res[i]*=ainv;
    }
    return res;
}

template<typename mint>
vector<mint> fps_inv(vector<mint> a){
    int n=(int)a.size();
    assert(n>0&&a[0]!=0);
    int cnt=nonzero_count(a);
    int thres=mint::can_ntt()?123:727;
    if(cnt<=thres) return fps_inv_sparse<mint>(a);
    return fps_inv_dense<mint>(a);
}
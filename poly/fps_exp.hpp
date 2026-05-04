#pragma once

#include "library/poly/convolution.hpp"
#include "library/mod/modint_basic.hpp"
#include "library/poly/fps_basic.hpp"
#include "library/poly/fps_log.hpp"

template<typename mint>
vector<mint> fps_exp_dense(vector<mint> a){
    if(a.empty()) return {};
    int n=a.size();
    assert(a[0]==mint(0));
    vector<mint> q(1,1),b=a;
    b[0]+=1;
    for(int m=1; m<n; m<<=1){
        if(n<m*2) b.resize(m*2);
        vector<mint> g=b,h=q;
        g.resize(m*2);
        h.resize(m*2),h=fps_log<mint>(h);
        for(int i=0; i<m*2; ++i) g[i]-=h[i];
        q=convolution<mint>(q,g);
        q.resize(m*2);
    }
    q.resize(n);
    return q;
}

template<typename mint>
vector<mint> fps_exp_sparse(vector<mint> a){
    if(a.empty()) return {};
    int n=(int)a.size();
    assert(a[0]==mint(0));
    vector<pair<int,mint>> b; // non-zero
    for(int i=1; i<n; ++i){
        if(a[i]!=mint(0)) b.push_back({i-1,a[i]*mint(i)});
    }
    vector<mint> res(n);
    res[0]=1;
    for(int i=1; i<n; ++i){
        mint val=0;
        for(auto [j,x]: b){
            if(j>i-1) break;
            val+=res[i-1-j]*x;
        }
        res[i]=val*inv<mint>(i);
    }
    return res;
}

template<typename mint>
vector<mint> fps_exp(vector<mint> a){
    if(a.empty()) return {};
    assert(a[0]==mint(0));
    int cnt=nonzero_count(a);
    int thres=mint::can_ntt()?234:2345;
    if(cnt<=thres) return fps_exp_sparse<mint>(a);
    return fps_exp_dense<mint>(a);
}
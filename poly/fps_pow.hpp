#pragma once

#include "library/poly/convolution.hpp"
#include "library/mod/modint_basic.hpp"
#include "library/poly/fps_basic.hpp"
#include "library/poly/fps_log.hpp"
#include "library/poly/fps_exp.hpp"

// const term = 1: fps_pow_1

template<typename mint>
vector<mint> fps_pow_1_dense(vector<mint> a, mint k){
    assert(!a.empty()&&a[0]==mint(1));
    a=fps_log<mint>(a);
    for(int i=0; i<(int)a.size(); ++i) a[i]*=k;
    a=fps_exp<mint>(a);
    return a;
}

template<typename mint>
vector<mint> fps_pow_1_sparse(vector<mint> a, mint k){
    assert(!a.empty()&&a[0]==mint(1));
    int n=(int)a.size();
    vector<pair<int,mint>> b; // non-zero
    for(int i=1; i<n; ++i){
        if(a[i]!=mint(0)) b.push_back({i,a[i]});
    }
    vector<mint> res(n);
    res[0]=1;
    for(int i=1; i<n; ++i){
        mint val=0;
        for(auto [j,x]: b){
            if(j>i) break;
            mint tmp=res[i-j]*x;
            val+=tmp*(mint(j)*k-mint(i-j));
        }
        res[i]=val*inv<mint>(i);
    }
    return res;
}

template<typename mint>
vector<mint> fps_pow_1(vector<mint> a, mint k){
    assert(!a.empty()&&a[0]==mint(1));
    int cnt=nonzero_count(a);
    int thres=mint::can_ntt()?123:1234;
    if(cnt<=thres) return fps_pow_1_sparse<mint>(a,k);
    return fps_pow_1_dense<mint>(a,k);
}

template<typename mint>
vector<mint> fps_pow(vector<mint> a, ll k){
    assert(k>=0);
    if(a.empty()){
        if(k==0) return {1};
        return {};
    }
    int n=a.size();
    if(k==0){
        vector<mint> res(n);
        res[0]=1;
        return res;
    }
    int m=n;
    for(int i=n-1; i>=0; --i){
        if(a[i]!=mint(0)) m=i;
    }
    if(m>0&&(k>=n||1ll*m*k>=n)){
        return vector<mint>(n);
    }
    int offset=m*k;
    mint c=a[m],cinv=c.inv();
    vector<mint> b(n-offset);
    for(int i=0; i<n-offset; ++i) b[i]=a[m+i]*cinv;
    b=fps_pow_1<mint>(b,mint(k));
    vector<mint> res(n);
    c=c.pow(k);
    for(int i=0; i<n-offset; ++i) res[offset+i]=b[i]*c;
    return res;
}
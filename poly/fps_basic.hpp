#pragma once

#include "library/mod/modint_basic.hpp"

template<typename mint>
vector<mint> derivative(vector<mint> a){
    if(a.empty()) return {};
    int n=a.size();
    vector<mint> res(n-1);
    for(int i=0; i<n-1; ++i) res[i]=a[i+1]*mint(i+1);
    return res;
}

template<typename mint>
vector<mint> integral(vector<mint> a){
    int n=a.size();
    vector<mint> res(n+1);
    for(int i=0; i<n; ++i) res[i+1]=a[i]*inv<mint>(i+1);
    return res;
}

template<typename mint>
vector<mint> dot(vector<mint> a, vector<mint> b){
    int n=a.size();
    assert(n==(int)b.size());
    vector<mint> res(n);
    for(int i=0; i<n; ++i) res[i]=a[i]*b[i];
    return res;
}

template<typename mint>
int nonzero_count(vector<mint> a){
    int res=0;
    for(int i=0; i<(int)a.size(); ++i) res+=(a[i]!=mint(0));
    return res;
}
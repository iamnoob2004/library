#pragma once

#include "library/poly/fps_inv.hpp"

template<typename mint>
pair<vector<mint>,vector<mint>> poly_div(vector<mint> a, vector<mint> b){
    assert(b.back()!=mint(0));
    int n=a.size(),m=b.size();
    if(n<m) return {{},a};
    int k=n-m+1;
    auto ra=a,rb=b;
    reverse(ra.begin(),ra.end());
    reverse(rb.begin(),rb.end());
    ra.resize(k),rb.resize(k);
    auto q=convolution<mint>(ra,fps_inv<mint>(rb));
    q.resize(k);
    reverse(q.begin(),q.end());
    auto c=convolution<mint>(q,b);
    for(int i=0; i<n; ++i) a[i]-=c[i];
    while(!a.empty()&&a.back()==mint(0)) a.pop_back();
    return {q,a};
}
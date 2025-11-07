#pragma once

#include "library/poly/ntt.hpp"

template<typename mint>
vector<mint> convolution_naive(vector<mint> a, vector<mint> b){
    if(a.empty()||b.empty()) return {};
    int n=((int)a.size())+((int)b.size())-1;
    vector<mint> res(n);
    for(int i=0; i<((int)a.size()); ++i) for(int j=0; j<((int)b.size()); ++j){
        res[i+j]+=a[i]*b[j];
    }
    return res;
}

template<typename mint>
vector<mint> convolution_ntt(vector<mint> a, vector<mint> b){
    if(a.empty()||b.empty()) return {};
    static NTT<mint> ntt;
    int n=((int)a.size())+((int)b.size())-1;
    int m=1,k=0;
    while(m<n) m<<=1,k++;
    a.resize(m),b.resize(m);
    ntt.trans(a,k),ntt.trans(b,k);
    for(int i=0; i<m; ++i) a[i]*=b[i];
    ntt.trans(a,k,true);
    a.resize(n);
    return a;
}

template<typename mint>
vector<mint> convolution(vector<mint> a, vector<mint> b){
    if(a.empty()||b.empty()) return {};
    int n=((int)a.size())+((int)b.size())-1;
    if(mint::ntt_data().first<0||n<49) return convolution_naive(a,b);
    return convolution_ntt(a,b);
}
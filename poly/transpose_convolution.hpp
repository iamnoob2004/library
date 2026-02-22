#pragma once

#include "library/poly/ntt.hpp"

template<typename mint>
vector<mint> transpose_convolution_naive(vector<mint> a, vector<mint> b){
    if(a.empty()||b.empty()) return {};
    int n=((int)a.size())-((int)b.size())+1;
    vector<mint> res(n);
    for(int i=0; i<n; ++i) for(int j=0; j<((int)b.size()); ++j){
        res[i]+=a[i+j]*b[j];
    }
    return res;
}

template<typename mint>
vector<mint> transpose_convolution_ntt(vector<mint> a, vector<mint> b){
    if(a.empty()||b.empty()) return {};
    static NTT<mint> ntt;
    int n=((int)a.size()),m=((int)b.size());
    int l=1,k=0;
    while(l<n) l<<=1,k++;
    a.resize(l),b.resize(l);
    ntt.transpose_dft(a,k,true),ntt.dft(b,k);
    for(int i=0; i<l; ++i) a[i]*=b[i];
    ntt.transpose_dft(a,k);
    a.resize(n-m+1);
    return a;
}

template<typename mint>
vector<mint> transpose_convolution(vector<mint> a, vector<mint> b){
    if(a.empty()||b.empty()) return {};
    int n=((int)a.size())-((int)b.size())+1;
    if(mint::ntt_data().first<0||n<=5||max((int)a.size(),(int)b.size())<49) return transpose_convolution_naive(a,b);
    return transpose_convolution_ntt(a,b);
}
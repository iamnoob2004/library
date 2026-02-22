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
    if(n>10&&n-2<=m/2){
        mint a_last=a.back(),b_last=b.back();
        a.pop_back(),b.pop_back();
        vector<mint> res=convolution_ntt(a,b);
        res.resize(n);
        res[n-1]+=a_last*b_last;
        for(int i=0; i<((int)a.size()); ++i) res[i+((int)b.size())]+=a[i]*b_last;
        for(int i=0; i<((int)b.size()); ++i) res[i+((int)a.size())]+=b[i]*a_last;
        return res;
    }
    a.resize(m),b.resize(m);
    ntt.dft(a,k),ntt.dft(b,k);
    for(int i=0; i<m; ++i) a[i]*=b[i];
    ntt.dft(a,k,true);
    a.resize(n);
    return a;
}

template<typename mint>
vector<mint> convolution(vector<mint> a, vector<mint> b){
    if(a.empty()||b.empty()) return {};
    int n=((int)a.size())+((int)b.size())-1;
    if(mint::ntt_data().first<0||n<49||min((int)a.size(),(int)b.size())<=5) return convolution_naive(a,b);
    return convolution_ntt(a,b);
}
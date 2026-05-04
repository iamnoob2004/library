#pragma once

#include "library/mod/modint.hpp"
#include "library/poly/ntt.hpp"
#include "library/mod/garner.hpp"

// n,m = a.size(),b.size() (n >= m)
// res_i = \sum_j a_{j+i}b_j

template<typename mint>
vector<mint> middle_product_naive(vector<mint> a, vector<mint> b){
    if(a.size()<b.size()||a.empty()||b.empty()) return {};
    int n=((int)a.size())-((int)b.size())+1;
    vector<mint> res(n);
    for(int i=0; i<n; ++i) for(int j=0; j<((int)b.size()); ++j){
        res[i]+=a[i+j]*b[j];
    }
    return res;
}

template<typename mint>
vector<mint> middle_product_ntt(vector<mint> a, vector<mint> b){
    if(mint::ntt_data().first<0) return {};
    if(a.size()<b.size()||a.empty()||b.empty()) return {};
    static NTT<mint> ntt;
    int n=((int)a.size())-((int)b.size())+1;
    int start=((int)b.size())-1;
    int m=1,k=0;
    while(m<((int)a.size())) m<<=1,k++;
    reverse(b.begin(),b.end());
    a.resize(m),b.resize(m);
    ntt.dft(a,k),ntt.dft(b,k);
    for(int i=0; i<m; ++i) a[i]*=b[i];
    ntt.dft(a,k,true);
    a=vector<mint>(a.begin()+start,a.begin()+start+n);
    return a;
}

template<typename mint>
vector<mint> middle_product_garner(vector<mint> a, vector<mint> b){
    if(a.size()<b.size()||a.empty()||b.empty()) return {};
    static constexpr int p0=167772161,p1=469762049,p2=754974721;
    using mint0=modint<p0>;
    using mint1=modint<p1>;
    using mint2=modint<p2>;
    int n=a.size(),m=b.size();
    vector<mint0> a0(n),b0(m);
    vector<mint1> a1(n),b1(m);
    vector<mint2> a2(n),b2(m);
    for(int i=0; i<n; ++i) a0[i]=a[i].x,a1[i]=a[i].x,a2[i]=a[i].x;
    for(int i=0; i<m; ++i) b0[i]=b[i].x,b1[i]=b[i].x,b2[i]=b[i].x;
    vector<mint0> res0=middle_product_ntt<mint0>(a0,b0);
    vector<mint1> res1=middle_product_ntt<mint1>(a1,b1);
    vector<mint2> res2=middle_product_ntt<mint2>(a2,b2);
    vector<mint> res(res0.size());
    for(int i=0; i<((int)res.size()); ++i){
        res[i]=garner3<mint>({res0[i].x,res1[i].x,res2[i].x},{p0,p1,p2});
    }
    return res;
}

template<typename mint>
vector<mint> middle_product(vector<mint> a, vector<mint> b){
    if(a.size()<b.size()||a.empty()||b.empty()) return {};
    int n=((int)a.size())-((int)b.size())+1;
    if(min(n,(int)b.size())<=5) return middle_product_naive(a,b);
    if(mint::ntt_data().first>=0){
        if(((int)a.size())<49) return middle_product_naive(a,b);
        return middle_product_ntt(a,b);
    }
    if(((int)a.size())<333) return middle_product_naive(a,b);
    return middle_product_garner(a,b);
}
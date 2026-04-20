#pragma once

#include "library/poly/convolution.hpp"

// calculate f(ar^k) for 0 <= k < m
template<typename mint>
vector<mint> multipoint_eval_on_geom_seq(vector<mint> f, mint a, mint r, int m){
    int n=f.size();
    if(m==0) return {};
    auto eval=[&](mint x) -> mint{
        mint res,pw=1;
        for(int i=0; i<n; ++i){
            res+=f[i]*pw;
            pw*=x;
        }
        return res;
    };
    if(r==mint(0)){
    	vector<mint> res(m);
    	for(int i=1; i<m; ++i) res[i]=f[0];
    	res[0]=eval(a);
    	return res;
    }
    if(n<49||m<49){
    	vector<mint> res(m);
    	for(int i=0; i<m; ++i,a*=r) res[i]=eval(a);
    	return res;
    }
    {
    	mint pw=1;
    	for(int i=0; i<n; ++i,pw*=a) f[i]*=pw;
    }
    auto calc=[&](mint q, int k){
        // calc q^binom(i,2) for 0 <= i < k
        vector<mint> res(k);
        mint pw=1;
        res[0]=1;
        for(int i=0; i<k-1; ++i,pw*=q){
            res[i+1]=res[i]*pw;
        }
        return res;
    };
    vector<mint> vec1=calc(r.inv(),max(n,m)),vec2=calc(r,n+m);
    vector<mint> c=f,d=vec2;
    for(int i=0; i<n; ++i){
        c[i]=f[n-1-i]*vec1[n-1-i];
    }
    c=convolution<mint>(c,vec2);
    vector<mint> res(m);
    for(int i=0; i<m; ++i){
        res[i]=c[n-1+i]*vec1[i];
    }
    return res;
}
#pragma once

#include "library/mod/modint_basic.hpp"
#include "library/poly/convolution.hpp"

// sum_j w_j[x^j]f^i for i=0,1,...,m
template<typename mint>
vector<mint> power_projection(vector<mint> f, vector<mint> w, int m){
    assert(f.size()==w.size());
    if(f.empty()) return vector<mint>(m+1,mint(0));
    // make [x^0]f = 0 with 1 convolution, faster than poly inv
    if(f[0]!=mint(0)){
        mint c=f[0];
        f[0]=0;
        vector<mint> res=power_projection(f,w,m);
        for(int i=0; i<=m; ++i) res[i]*=ifac<mint>(i);
        vector<mint> tmp(m+1);
        mint pow_c=1;
        for(int i=0; i<=m; ++i) tmp[i]=pow_c*ifac<mint>(i),pow_c*=c;
        res=convolution<mint>(res,tmp);
        res.resize(m+1);
        for(int i=0; i<=m; ++i) res[i]*=fac<mint>(i);
        return res;
    }
    int n=1,k=1;
    while(n<(int)f.size()) n<<=1;
    f.resize(n),w.resize(n);
    reverse(w.begin(),w.end());
    // store in y-major order
    // omit y^k for Q
    vector<mint> P(2*n),Q(2*n); 
    for(int i=0; i<n; ++i) P[i]=w[i],Q[i]=-f[i];
    while(n>1){
        vector<mint> R(2*n*k); // Q(-x)
        for(int i=0; i<2*n*k; ++i) R[i]=(i%2==0?Q[i]:-Q[i]);
        // P(x,y)Q(-x,y),Q(x,y)Q(-x,y)
        vector<mint> PQ=convolution(P,R),QQ=convolution(Q,R);
        PQ.resize(4*n*k),QQ.resize(4*n*k);
        // y^k
        for(int i=0; i<2*n*k; ++i){
            PQ[2*n*k+i]+=P[i];
            QQ[2*n*k+i]+=Q[i]+R[i];
        }
        fill(P.begin(),P.end(),mint(0)),fill(Q.begin(),Q.end(),mint(0));
        for(int j=0; j<2*k; ++j) for(int i=0; i<n/2; ++i){
            // x^i*y^j
            P[n*j+i]=PQ[2*n*j+2*i+1];
            Q[n*j+i]=QQ[2*n*j+2*i];
        }
        n/=2,k*=2;
    }
    // reverse y
    vector<mint> p(k);
    for(int i=0; i<k; ++i) p[i]=P[2*i];
    reverse(p.begin(),p.end());
    p.resize(m+1);
    return p;
}
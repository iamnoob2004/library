#pragma once

#include "library/poly/poly.hpp"
#include "library/poly/transpose_convolution.hpp"

template<typename mint>
poly<mint> composition(poly<mint> f, poly<mint> g){
    int _n=max<int>(f.size(),g.size());
    int n=1;
    while(n<_n) n<<=1;
    f.resize(n),g.resize(n);
    auto dfs=[&](auto &self, int n, int k, vector<mint> Q) -> vector<mint>{
        if(n==1){
            vector<mint> p(2*k);
            reverse(f.begin(),f.end());
            for(int i=0; i<k; ++i) p[2*i]=f[i];
            return p;
        }
        vector<mint> R(2*n*k); // Q(-x)
        for(int i=0; i<2*n*k; ++i) R[i]=(i%2==0?Q[i]:-Q[i]);
        vector<mint> QQ=convolution(Q,R);
        QQ.resize(4*n*k);
        // y^k
        for(int i=0; i<2*n*k; ++i) QQ[2*n*k+i]+=Q[i]+R[i];
        vector<mint> nxt_Q(2*n*k);
        for(int j=0; j<2*k; ++j) for(int i=0; i<n/2; ++i){
            // x^i*y^j
            nxt_Q[n*j+i]=QQ[2*n*j+2*i];
        }

        vector<mint> nxt_p=self(self,n/2,k*2,nxt_Q);
        vector<mint> pq(4*n*k);
        for(int j=0; j<2*k; ++j) for(int i=0; i<n/2; ++i){
            // x^i*y^j
            pq[2*n*j+2*i+1]=nxt_p[n*j+i];
        }
        vector<mint> p(2*n*k);
        for(int i=0; i<2*n*k; ++i) p[i]=pq[2*n*k+i];
        pq.pop_back();
        vector<mint> x=transpose_convolution(pq,R);
        for(int i=0; i<2*n*k; ++i) p[i]+=x[i];
        return p;
    };
    int k=1;
    vector<mint> Q(2*n);
    for(int i=0; i<n; ++i) Q[i]=-g[i];
    vector<mint> p=dfs(dfs,n,k,Q);
    vector<mint> res(n);
    for(int i=0; i<n; ++i) res[i]=p[i];
    reverse(all(res));
    res.resize(_n);
    return res;
}
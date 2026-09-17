#pragma once

#include "poly/transpose_convolution.hpp"

template<typename mint>
vector<mint> composition_small(vector<mint> f, vector<mint> g){
    int n=max<int>(f.size(),g.size());
    if(f.empty()){
        return vector<mint>(n,0);
    }
    int k=2;
    while(k*k<n) k++;
    vector<vector<mint>> vec1(k+1),vec2(k);
    vec1[0]={1};
    for(int i=1; i<=k; ++i){
        vec1[i]=convolution<mint>(vec1[i-1],g);
        if((int)vec1[i].size()>n) vec1[i].resize(n);
    }
    vec2[0]={1},vec2[1]=vec1[k];
    for(int i=2; i<k; ++i){
        vec2[i]=convolution<mint>(vec2[i-1],vec2[1]);
        if((int)vec2[i].size()>n) vec2[i].resize(n);
    }
    vector<vector<mint>> vec3(k);
    for(int i=0; i<k; ++i){
        vec3[i].resize(n);
        for(int j=0; j<k; ++j){
            if(i*k+j>=(int)f.size()) break;
            mint w=f[i*k+j];
            for(int t=0; t<(int)vec1[j].size(); ++t){
                vec3[i][t]+=w*vec1[j][t];
            }
        }
    }
    vector<mint> res(n);
    for(int i=0; i<k; ++i){
        vector<mint> tmp=convolution<mint>(vec3[i],vec2[i]);
        tmp.resize(n);
        for(int j=0; j<n; ++j) res[j]+=tmp[j];
    }
    return res;
}

template<typename mint>
vector<mint> composition_big(vector<mint> f, vector<mint> g){
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

template<typename mint>
vector<mint> composition(vector<mint> f, vector<mint> g){
    int n=max<int>(f.size(),g.size());
    if(n<=1000) return composition_small<mint>(f,g);
    return composition_big<mint>(f,g);
}
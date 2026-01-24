#pragma once

#include "library/poly/poly.hpp"

template<typename T>
T bostan_mori(poly<T> P, poly<T> Q, ll k){
    // [x^k] a/b
    for(; k; k>>=1){
        poly<T> Q_neg=Q;
        for(int i=0; i<(int)Q_neg.size(); ++i){
            if(i&1) Q_neg[i]=-Q_neg[i];
        }
        P*=Q_neg;
        Q*=Q_neg;
        for(int i=1; i<(int)Q.size(); i+=2){
            assert(Q[i]==T(0));
        }
        poly<T> nw_P,nw_Q;
        for(int i=0; i<(int)Q.size(); i+=2){
            nw_Q.push_back(Q[i]);
        }
        for(int i=k&1; i<(int)P.size(); i+=2){
            nw_P.push_back(P[i]);
        }
        P=nw_P,Q=nw_Q;
    }
    return P[0]/Q[0];
}

template<typename T>
pair<poly<T>,poly<T>> get_genfunc(vector<T> a, vector<T> b){
    // linear recursion to rational function
    // a: base case
    // a_m = \sum_{i=0}^{n-1} b_i*a_{m-1-i}
    int n=a.size();
    assert(a.size()==b.size());
    poly<T> Pa(a),Pb(b);
    Pb.insert(Pb.begin(),T(0));
    poly<T> P=Pa*Pb-Pa,Q=Pb;
    P.resize(n);
    Q[0]-=1;
    return {P,Q};
}

template<typename T>
T linear_recursion_kth(vector<T> a, vector<T> b, ll k){
    int n=a.size();
    assert(a.size()==b.size());
    if(k<0) return 0;
    if(k<n) return a[k];
    auto [P,Q]=get_genfunc(a,b);
    // a = P/Q
    return bostan_mori(P,Q,k);
}

/*
template<typename T>
vector<T> linear_recursion_consecutive(vector<T> a, vector<T> b, ll k, int m){
    // a: base case
    // a_m = \sum_{i=0}^{n-1} b_i*a_{m-1-i}
    // return k-th to (k+m-1)-th term
    int n=a.size();
    assert(a.size()==b.size());
    poly<T> Pa(a),Pb(b);
    Pb.insert(Pb.begin(),T(0));
    poly<T> diff=Pa*Pb-Pa;
    diff.resize(n);
    Pb[0]-=1;
    Pb.resize(n+m);
    diff*=Pb.inverse();
    a=vector<T>(diff);
    if(k<n){
        vector<T> ans(m);
        for(int i=0; i<m; ++i){
            if(i+k>=0) ans[i]=a[i+k];
        }
        return ans;
    }
    poly<T> tar(n+1,1);
    for(int i=0; i<n; ++i) tar[i]=-b[n-1-i];
    auto get=[&](auto &self, ll _k) -> poly<T>{
        if(_k==0) return poly<T>({1});
        poly<T> P=self(self,_k>>1);
        P*=P;
        P=P.divide(tar).second;
        if(_k&1){
            P*=poly<T>({0,1});
            P=P.divide(tar).second;
        }
        return P;
    };
    poly<T> res=get(get,k);
    assert((int)res.size()==n);
    reverse(res.begin(),res.end());
    res*=diff;
    vector<T> ans(m);
    for(int i=0; i<m; ++i) ans[i]=res[i+n-1];
    return ans;
}
*/
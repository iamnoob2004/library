#pragma once

#include "poly/fps.hpp"

template<typename T, typename S>
T bostan_mori(fps<T> P, fps<T> Q, S k){
    // [x^k] a/b
    for(; k; k>>=1){
        fps<T> Q_neg=Q;
        for(int i=0; i<(int)Q_neg.size(); ++i){
            if(i&1) Q_neg[i]=-Q_neg[i];
        }
        P*=Q_neg;
        Q*=Q_neg;
        for(int i=1; i<(int)Q.size(); i+=2){
            assert(Q[i]==T(0));
        }
        fps<T> nw_P,nw_Q;
        for(int i=0; i<(int)Q.size(); i+=2){
            nw_Q.push_back(Q[i]);
        }
        for(int i=int(k&S(1)); i<(int)P.size(); i+=2){
            nw_P.push_back(P[i]);
        }
        P=nw_P,Q=nw_Q;
    }
    return P[0]/Q[0];
}

template<typename T>
pair<fps<T>,fps<T>> get_genfunc(vector<T> a, vector<T> b){
    // linear recursion to rational function
    // a: base case
    // a_m = \sum_{i=0}^{n-1} b_i*a_{m-1-i}
    int n=a.size();
    assert(a.size()==b.size());
    fps<T> Pa(a),Pb(b);
    Pb.insert(Pb.begin(),T(0));
    fps<T> P=Pa*Pb-Pa,Q=Pb;
    P.resize(n);
    Q[0]-=1;
    return {P,Q};
}

template<typename T, typename S>
T linear_recursion_kth(vector<T> a, vector<T> b, S k){
    int n=a.size();
    assert(a.size()==b.size());
    if(k<S(0)) return 0;
    if(k<S(n)) return a[(int)k];
    auto [P,Q]=get_genfunc(a,b);
    // a = P/Q
    return bostan_mori(P,Q,k);
}
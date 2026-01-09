#pragma once

#include "library/linalg/arbitrary_matrix.hpp"

template<typename T>
matrix<T> solve_equation(vector<matrix<mint>> a){
    int n=a.size();
    for(int i=0; i<n; ++i){
        assert(a[i].h==1&&a[i].w==n+1);
    }
    for(int i=0; i<n; ++i){
        if(a[i][i]==T(0)){
            for(int j=i+1; j<n; ++j) if(a[j][i]!=T(0)){
                swap(a[i],a[j]);
                break;
            }
            if(a[i][i]==T(0)) return matrix<T>();
        }
        T x=T(1)/a[i][i];
        a[i]*=x;
        assert(a[i][i]==T(1));
        for(int j=0; j<n; ++j) if(i!=j) {
            x=-a[j][i];
            a[j]+=a[i]*x;
        }
    }
    matrix<T> res;
    res.set_hw(n,1);
    for(int i=0; i<n; ++i) res[i]=-a[i][n];
    return res;
}
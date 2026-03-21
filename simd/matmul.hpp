#pragma once

#include "library/simd/dot_product_mod.hpp"

#ifndef i_am_noob

// b transposed
vector<int> matmul_mod(int n, int m, int k, vector<int> a, vector<int> b, int mod){
    assert((int)a.size()==n*k&&(int)b.size()==m*k);
    vector<int> res(n*m);
    int kk=k+1;
    while(kk%4) kk++;
    kk>>=2;
    __m256i vec_a[n*kk],vec_b[m*kk];
    for(int i=0; i<n; ++i){
        ll tmp[kk*4]{};
        for(int j=0; j<k; ++j){
            tmp[j]=a[i*k+j];
        }
        for(int j=0; j<kk; ++j){
            vec_a[i*kk+j]=_mm256_loadu_si256((__m256i*)&tmp[j*4]);
        }
    }
    for(int i=0; i<m; ++i){
        ll tmp[kk*4]{};
        for(int j=0; j<k; ++j){
            tmp[j]=b[i*k+j];
        }
        for(int j=0; j<kk; ++j){
            vec_b[i*kk+j]=_mm256_loadu_si256((__m256i*)&tmp[j*4]);
        }
    }
    for(int i=0; i<n; ++i){
        for(int j=0; j<m; ++j){
            res[i*m+j]=dot_product_mod(kk,vec_a+(i*kk),vec_b+(j*kk),mod);
        }
    }
    return res;
}

#else

vector<int> matmul_mod(int n, int m, int k, vector<int> a, vector<int> b, int mod){
    return {};
}

#endif
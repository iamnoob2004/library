#pragma once

#include "library/ds/dynamic_bitset.hpp"

int bitset_lcs(vector<int> a, vector<int> b, int mx){
    mx++;
    int n=a.size();
    vector<dynamic_bitset> vec(mx,dynamic_bitset(n+1));
    for(int i=0; i<n; ++i) vec[a[i]].set(i);

    dynamic_bitset res(n+1),tmp(n+1);
    for(int y: b){
        tmp=vec[y];
        tmp|=res;
        res<<=1;
        res.set(0);
        res=tmp-res;
        res.flip();
        res&=tmp;
    }
    return res.count();
}
#pragma once

template<typename mint>
vector<mint> pow_table(mint a, int n){
    vector<mint> res(n,1);
    for(int i=1; i<n; ++i) res[i]=res[i-1]*a;
    return res;
}
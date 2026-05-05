#pragma once

template<typename mint>
vector<mint> all_inverse(vector<mint> a){
    int n=a.size();
    for(int i=0; i<n; ++i) assert(a[i]!=mint(0));
    vector<mint> pref(n+1);
    pref[0]=1;
    for(int i=0; i<n; ++i) pref[i+1]=pref[i]*a[i];
    mint cur=pref[n].inv();
    vector<mint> res(n);
    for(int i=n-1; i>=0; --i){
        res[i]=pref[i]*cur;
        cur*=a[i];
    }
    return res;
}
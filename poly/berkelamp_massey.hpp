#pragma once

// s = {1,1,2,3,5,8,13} -> returns {1,1}
template<typename T>
vector<T> berlekamp_massey(const vector<T> &s) {
    int n=s.size(),l=0,m=1;
    vector<T> b(n),c(n);
    T ld=b[0]=c[0]=T(1);
    for(int i=0; i<n; i++,m++){
        T d=s[i];
        for(int j=1; j<=l; j++) d+=c[j]*s[i-j];
        if(d==0) continue;
        vector<T> tmp=c;
        T coef=d/ld;
        for(int j=m; j<n; j++) c[j]-=coef*b[j-m];
        if(2*l<=i){
            l=i+1-l,ld=d,m=0;
            b=tmp;
        }
    }
    c.resize(l+1);
    c.erase(c.begin());
    for(T &x: c) x=-x;
    return c;
}
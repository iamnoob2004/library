#pragma once

int get_mex(vector<int> a){
    int n=a.size();
    vector<int> res(n+1);
    for(auto i: a){
        if(i<=n) res[i]=1;
    }
    for(int i=0; i<=n; ++i){
        if(!res[i]) return i;
    }
    return -1;
}
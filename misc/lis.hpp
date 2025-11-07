#pragma once

template<typename T>
int LIS(vector<T> a){
    vector<T> vec;
    for(T x: a){
        auto it=lower_bound(vec.begin(),vec.end(),x);
        if(it==vec.end()) vec.pb(x);
        else (*it)=x;
    }
    return vec.size();
}
#pragma once

template<typename T>
struct chtholly_tree{
    ll n;
    map<ll,T> mp;

    chtholly_tree():chtholly_tree(0){}
    chtholly_tree(ll _n){
        n=_n;
        mp[0]=mp[n]=T();
    }
    
    // returns overwritten data
    vector<pair<pair<ll,ll>,T>> assign(ll l, ll r, T x){
        assert(0<=l&&l<=r&&r<=n);
        if(l==r) return {};
        auto it=mp.lower_bound(l);
        if(it->first>l){
            mp[l]=prev(it)->second;
        }
        it=mp.lower_bound(r);
        if(it->first>r){
            mp[r]=prev(it)->second;
        }
        it=mp.lower_bound(l);
        vector<pair<pair<ll,ll>,T>> res;
        while(it->first<r){
            res.pb({{it->first,next(it)->first},it->second});
            it++;
            mp.erase(prev(it));
        }
        mp[l]=x;
        return res;
    }
};
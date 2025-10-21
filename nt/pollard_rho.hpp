#pragma once

#include "library/nt/binary_gcd.hpp"
#include "library/nt/miller_rabin.hpp"

// https://judge.yosupo.jp/problem/factorize
struct pollard_rho{
    ll find_factor(ll n){
        if(n<=1||is_prime(n)) return 1;
        static const vector<int> small={2,3,5,7,11,13,17,19};
        for(int p: small){
            if(n%p==0) return p;
        }
        using mint=montgomery_modint_64<777771449,false>;
        mint::set_mod(n);
        mint x,y(2),d,t(1);
        auto f=[&](mint a){return a*a+t;};
        for(int l=2; ; l<<=1){
            x=y;
            int m=min(l,32);
            for(int i=0; i<l; i+=m){
                d=1;
                for(int j=0; j<m; ++j){
                    y=f(y);
                    d*=x-y;
                }
                ll g=bgcd<ll>(d.get(),n);
                if(g==n){
                    l=1,y=2,t+=1;
                    break;
                }
                if(g!=1) return g;
            }
        }
    }
    map<ll,int> mp;
    void dfs(ll n){
        if(n<=1) return;
        if(is_prime(n)) return mp[n]++,void();
        ll d=find_factor(n);
        dfs(d);
        dfs(n/d);
    }
};

vector<pair<ll,int>> factorize(ll n){
    pollard_rho tmp;
    tmp.dfs(n);
    vector<pair<ll,int>> res;
    for(auto [x,y]: tmp.mp){
        res.pb({x,y});
    }
    return res;
}

vector<ll> find_all_divisors(ll n){
    vector<pair<ll,int>> vec=factorize(n);
    vector<ll> res;
    auto dfs=[&](auto &self, int i, ll cur){
        if(i==(int)vec.size()){
            res.pb(cur);
            return;
        }
        for(int j=0; j<vec[i].second; ++j){
            self(self,i+1,cur);
            cur*=vec[i].first;
        }
        self(self,i+1,cur);
    };
    dfs(dfs,0,1);
    sort(res.begin(),res.end());
    return res;
}
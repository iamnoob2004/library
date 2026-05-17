#pragma once

// n=10: 0,1,2,3,5,10
struct array_of_floor{
    ll n;
    int m,t;

    array_of_floor(ll _n=0){init(_n);}

    void init(ll _n){
        n=_n;
        t=0;
        while(1ll*(t+1)*(t+1)<=n) t++;
        m=t*2+1;
        if(1ll*t*(t+1)>n) m--;
    }

    int get_id(ll x){
        if(x<=t) return x;
        return m-n/x;
    }

    ll operator [](int i){
        if(i<=t) return i;
        return n/(m-i);
    }
};
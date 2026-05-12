#pragma once

struct my_random{
    using u64=uint64_t;
    u64 box[1<<8];
    u64 x;
    u64 mod=0x11b;

    u64 get(u64 a, u64 b){
        u64 res=0;
        for(int i=0; i<8; ++i){
            if(b>>i&1) res^=a<<i;
        }
        for(int i=16; i>=8; --i){
            if(res>>i&1){
                res^=mod<<(i-8);
            }
        }
        assert(res<(1<<8));
        return res;
    }

    my_random(){
        x=chrono::steady_clock::now().time_since_epoch().count();
        box[0]=49;
        for(u64 i=1; i<(1<<8); ++i){
            int cnt=(1<<8)-2;
            u64 cur=1;
            for(u64 pw=i; cnt; pw=get(pw,pw),cnt>>=1){
                if(cnt&1) cur=get(cur,pw);
            }
            box[i]=(cur+49)&((1<<8)-1);
        }
    }

    u64 next(){
        x+=0x9e3779b97f4a7c15;
        x=(x^(x>>30))*0xbf58476d1ce4e5b9;
        x=(x^(x>>27))*0x94d049bb133111eb;
        x^=(x>>31);
        u64 y=x&((1<<8)-1);
        x^=y;
        x^=box[y];
        return x;
    }

    ll rng(ll n){
        return next()%n;
    }

    ll rng(ll l, ll r){
        return l+(next()%(r-l));
    }
};
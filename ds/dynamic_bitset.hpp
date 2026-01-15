#pragma once

struct dynamic_bitset{
    using u64=uint64_t;

    int n,m;
    vector<u64> vec;

    void init(int _n){
        m=((_n+63)>>6);
        n=m<<6;
        vec.reserve(m);
        vec.assign(m,0);
    }
    dynamic_bitset(int _n=0){
        init(_n);
    }

    bool operator [](int i){return vec[i>>6]>>(i&63)&1;}
    void set(int i){vec[i>>6]|=u64(1)<<(i&63);}
    void reset(int i){vec[i>>6]&=~(u64(1)<<(i&63));}
    void flip(int i){vec[i>>6]^=u64(1)<<(i&63);}

    void flip(){
        for(int i=0; i<m; ++i) vec[i]=~vec[i];
    }

    dynamic_bitset operator &= (const dynamic_bitset &o){
        assert(m==o.m);
        for(int i=0; i<m; ++i) vec[i]&=o.vec[i];
        return *this;
    }
    dynamic_bitset operator |= (const dynamic_bitset &o){
        assert(m==o.m);
        for(int i=0; i<m; ++i) vec[i]|=o.vec[i];
        return *this;
    }
    dynamic_bitset operator ^= (const dynamic_bitset &o){
        assert(m==o.m);
        for(int i=0; i<m; ++i) vec[i]^=o.vec[i];
        return *this;
    }
    dynamic_bitset operator -= (const dynamic_bitset &o){
        assert(m==o.m);
        bool f=0,nf;
        for(int i=0; i<m; ++i){
            nf=o.vec[i]>vec[i]||(o.vec[i]==vec[i]&&f);
            if(f) vec[i]--;
            vec[i]-=o.vec[i];
            f=nf;
        }
        return *this;
    }

    dynamic_bitset operator - (const dynamic_bitset &o) const {return dynamic_bitset(*this)-=o;}

    dynamic_bitset operator <<= (int k){
        k=min(k,n);
        if(k>>6){
            int t=k>>6;
            for(int i=m-1; i>=t; --i) vec[i]=vec[i-t];
            for(int i=0; i<t; ++i) vec[i]=0;
            k-=t<<6;
        }
        if(k==0) return *this;
        for(int i=m-1; i>=0; --i){
            vec[i]<<=k;
            if(i) vec[i]|=vec[i-1]>>(64-k);
        }
        return *this;
    }

    dynamic_bitset operator >>= (int k){
        k=min(k,n);
        if(k>>6){
            int t=k>>6;
            for(int i=0; i<m-t; --i) vec[i]=vec[i+t];
            for(int i=m-t; i<m; ++i) vec[i]=0;
            k-=t<<6;
        }
        if(k==0) return *this;
        for(int i=0; i<m; ++i){
            vec[i]>>=k;
            if(i+1<m) vec[i]|=vec[i+1]<<(64-k);
        }
        return *this;
    }

    int find_next(int i){ // >= i
        u64 x=vec[i>>6]>>(i&63);
        if(x){
            return __builtin_ctzll(x)+i;
        }
        for(int j=(i>>6)+1; j<m; ++j) if(vec[j]){
            return __builtin_ctzll(vec[j])+(j<<6);
        }
        return n;
    }
    int find_prev(int i){ // <= i
        u64 x=vec[i>>6]<<(63-(i&63));
        if(x){
            return i-__builtin_clzll(x);
        }
        for(int j=(i>>6)-1; j>=0; --j) if(vec[j]){
            return (j<<6)+63-__builtin_clzll(vec[j]);
        }
        return -1;
    }
    int find_first(){
        return find_next(0);
    }
    int find_last(){
        return find_prev(n-1);
    }

    int count(){
        int res=0;
        for(int i=0; i<m; ++i) res+=__builtin_popcountll(vec[i]);
        return res;
    }
};
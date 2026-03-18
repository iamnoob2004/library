#pragma once

#include "library/nt/extgcd.hpp"

template<int id>
struct barrett_modint{
    using mint=barrett_modint;
    using u32=uint32_t;
    using u64=uint64_t;
    using u128=unsigned __int128;

    inline static u32 m;
    inline static u64 im;

    static u32 get(u64 a){
        // return a mod m
        u64 q=(u128(im)*a)>>64,r=a-q*m;
        if(r>=m) r-=m;
        return r;
    }

    static void set_mod(int _m){
        m=_m;
        im=((u64)(-1))/m;
    }
    static int get_mod(){
        return m;
    }

    u32 x;
    barrett_modint():x(0){}
    barrett_modint(int64_t _x):x((_x%m+m)%m){}

    mint operator += (const mint &o){
        x+=o.x;
        if(x>=m) x-=m;
        return *this;
    }
    mint operator -= (const mint &o){
        x+=m-o.x;
        if(x>=m) x-=m;
        return *this;
    }
    mint operator *= (const mint &o){
        x=mint::get(u64(x)*o.x);
        return *this;
    }
    mint operator /= (const mint &o){
        return (*this)*=o.inv();
    }
    mint operator + (const mint &o) const {return mint(*this)+=o;}
    mint operator - (const mint &o) const {return mint(*this)-=o;}
    mint operator * (const mint &o) const {return mint(*this)*=o;}
    mint operator / (const mint &o) const {return mint(*this)/=o;}
    mint operator - () const {return mint(0)-*this;}
    mint pow(int64_t n) const {
        assert(n>=0);
        mint res=1,b=*this;
        for(; n; n>>=1,b*=b) if(n&1) res*=b;
        return res;
    }
    inline mint inv1() const {
        return pow(m-2);
    }
    inline mint inv2() const {
        auto [g,val1,val2]=extgcd<int>(x,m);
        assert(g==1);
        return mint(val1);
    }
    mint inv() const {
        return inv2();
    }

    bool operator == (const mint &o) const {
        return x==o.x;
    }
    bool operator != (const mint &o) const {
        return x!=o.x;
    }

    friend istream& operator >> (istream& is, mint &b){
        int64_t y;
        is >> y;
        b=mint(y);
        return is;
    }
    friend ostream& operator << (ostream& os, const mint &b){
        return os << b.x;
    }

    // v2(m-1), 2^(v2(m-1))-th root
    static constexpr pair<int,int> ntt_data(){
        return {-1,-1};
    }
};

template<int id> using modint=barrett_modint<id>;

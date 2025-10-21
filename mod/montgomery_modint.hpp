#pragma once

#include "library/nt/extgcd.hpp"

// arbitrary modint, odd mod
// stores x*(2^K) mod m
// https://judge.yosupo.jp/problem/binomial_coefficient_prime_mod
// https://judge.yosupo.jp/problem/primality_test (used by miller rabin)
template<int id, bool is_prime, int K, typename word, typename dword, typename signed_word> // support multiple modulos at the same time
struct montgomery_modint{
    using mint=montgomery_modint;

    inline static word m,r,val64,m2; // m = modulo < 2^(K-2), r = (-m^(-1)) (mod 2^K), val64 = (2^(2K)) (mod m), m2 = 2m

    static void set_mod(word _m){
        assert((_m&1)&&_m<(word(1)<<(K-2)));
        m=_m,r=m,val64=(-dword(m))%m,m2=m*2;
        // use Newton's method to calculate p^(-1) (mod 2^K)
        // starts from p = p^(-1) (mod 4)
        for(int i=0; i<5; ++i) r*=2-m*r;
        r=-r;
        assert(r*m==word(-1));
    }
    static int get_mod(){
        return m;
    }

    word x;

    montgomery_modint():x(0){}
    montgomery_modint(int64_t _x):x(reduce(dword((_x%m+m)%m)*val64)){}

    word reduce(const dword &y) const {
        // (y + (yr mod 2^K)*p) / (2^K)
        // 0 <= return < 2p
        return (y+dword(word(y)*r)*m)>>K;
    }

    mint operator += (const mint &o){
        x+=o.x;
        if(x>=m2) x-=m2;
        return *this;
    }
    mint operator -= (const mint &o){
        x-=o.x;
        if(int32_t(x)<0) x+=m2;
        return *this;
    }
    mint operator *= (const mint &o){
        x=reduce(dword(x)*o.x);
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
        auto [g,val1,val2]=extgcd<signed_word>(get(),m);
        assert(g==1);
        return mint(val1);
    }
    mint inv() const {
        if(is_prime) return inv1();
        return inv2();
    }

    bool operator == (const mint &o) const {
        return (x>=m?x-m:x)==(o.x>=m?o.x-m:o.x);
    }
    bool operator != (const mint &o) const {
        return (x>=m?x-m:x)!=(o.x>=m?o.x-m:o.x);
    }

    word get() const {
        word res=reduce(x);
        return res>=m?res-m:res;
    }

    friend istream& operator >> (istream& is, mint &b){
        int64_t y;
        is >> y;
        b=mint(y);
        return is;
    }
    friend ostream& operator << (ostream& os, const mint &b){
        return os << b.get();
    }
};

template<int id, bool is_prime> using montgomery_modint_32=montgomery_modint<id,is_prime,32,uint32_t,uint64_t,int32_t>;
template<int id, bool is_prime> using montgomery_modint_64=montgomery_modint<id,is_prime,64,uint64_t,unsigned __int128,int64_t>;
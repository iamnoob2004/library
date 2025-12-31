#pragma once

#include "library/nt/extgcd.hpp"

template<int64_t m, bool is_prime, typename word, typename dword>
struct Modint{
    using mint=Modint;

    static word get_mod(){
        return m;
    }

    word x;
    Modint():x(0){}
    Modint(int64_t _x):x((_x%m+m)%m){}

    mint operator += (const mint &o){
        x+=o.x;
        if(x>=m) x-=m;
        return *this;
    }
    mint operator -= (const mint &o){
        x-=o.x;
        if(x<0) x+=m;
        return *this;
    }
    mint operator *= (const mint &o){
        x=((dword)x)*o.x%m;
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
        auto [g,val1,val2]=extgcd<word>(x,m);
        assert(g==1);
        return mint(val1);
    }
    mint inv() const {
        if(is_prime) return inv1();
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
        if(m==998244353) return {23,31};
        return {-1,-1};
    }
};

template<int64_t m, bool is_prime=true> using modint=Modint<m,is_prime,int32_t,int64_t>;
template<int64_t m, bool is_prime=true> using modint_64=Modint<m,is_prime,int64_t,__int128>;
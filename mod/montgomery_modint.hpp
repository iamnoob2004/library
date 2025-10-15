#pragma once

// arbitrary modint
// stores x*(2^32) mod p
// https://judge.yosupo.jp/problem/binomial_coefficient_prime_mod
template<int id> // support multiple modulos at the same time
struct montgomery_modint{
    using u32=uint32_t;
    using u64=uint64_t;
    using mint=montgomery_modint;

    static const int K=32;
    inline static u32 p,r,val64,p2; // p = modulo < 2^(K-2), r = (-p^(-1)) (mod 2^K), val64 = (2^(2K)) (mod p), p2 = 2p

    static void set_mod(u32 _p){
        assert((_p&1)&&_p<(u32(1)<<(K-2)));
        p=_p,r=p,val64=(-u64(p))%p,p2=p*2;
        // use Newton's method to calculate p^(-1) (mod 2^K)
        // starts from p = p^(-1) (mod 4)
        for(int i=0; i<4; ++i) r*=2-p*r;
        r=-r;
        assert(r*p==u32(-1));
    }
    static int get_mod(){
        return p;
    }

    u32 x;

    montgomery_modint():x(0){}
    montgomery_modint(int64_t _x):x(reduce(u64((_x%p+p)%p)*val64)){}

    u32 reduce(const u64 &y) const {
        // (y + (yr mod 2^K)*p) / (2^K)
        // 0 <= return < 2p
        return (y+u64(u32(y)*r)*p)>>K;
    }

    mint operator += (const mint &o){
        x+=o.x;
        if(x>=p2) x-=p2;
        return *this;
    }
    mint operator -= (const mint &o){
        x-=o.x;
        if(int32_t(x)<0) x+=p2;
        return *this;
    }
    mint operator *= (const mint &o){
        x=reduce(u64(x)*o.x);
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
    mint inv() const {
        return pow(p-2);
    }

    bool operator == (const mint &o) const {
        return (x>=p?x-p:x)==(o.x>=p?o.x-p:x);
    }
    bool operator != (const mint &o) const {
        return (x>=p?x-p:x)!=(o.x>=p?o.x-p:x);
    }

    u32 get() const {
        u32 res=reduce(x);
        return res>=p?res-p:res;
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


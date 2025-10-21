#pragma once

template<int m>
struct modint_composite{
    using mint=modint_composite;
    using u32=uint32_t;

    static int get_mod(){
        return m;
    }

    int x;
    modint_composite():x(0){}
    modint_composite(int64_t _x):x((_x%m+m)%m){}

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
        x=((int64_t)x)*o.x%m;
        return *this;
    }
    mint operator + (const mint &o) const {return mint(*this)+=o;}
    mint operator - (const mint &o) const {return mint(*this)-=o;}
    mint operator * (const mint &o) const {return mint(*this)*=o;}
    mint operator - () const {return mint(0)-*this;}
    mint pow(int64_t n) const {
        assert(n>=0);
        mint res=1,b=*this;
        for(; n; n>>=1,b*=b) if(n&1) res*=b;
        return res;
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
};
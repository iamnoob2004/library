#pragma once

template<int p>
struct modint{
    using u32=uint32_t;

    int x;
    modint():x(0){}
    modint(int64_t _x):x((_x%p+p)%p){}

    mint operator += (const mint &o){
        x+=o.x;
        if(x>=p) x-=p;
        return *this;
    }
    mint operator -= (const mint &o){
        x-=o.x;
        if(x<0) x+=p;
        return *this;
    }
    mint operator *= (const mint &o){
        x=((long long)x)*o.x%p;
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
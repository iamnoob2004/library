#pragma once

#include "library/nt/binary_gcd.hpp"

template<typename T>
struct frac{
    T abs(T x){
        if(x<0) x=-x;
        return x;
    }
    T p,q;
    void fix(){
        T g=bgcd<T>(abs(p),abs(q));
        if(q<0) g*=-1;
        p/=g,q/=g;
    }
    frac():p(0),q(1){}
    frac(T x):p(x),q(1){}
    frac(T x, T y):p(x),q(y){
        if(q<0) p=-p,q=-q;
    }
    frac operator +(const frac &o) const {return frac(p*o.q+q*o.p,q*o.q);}
    frac operator -(const frac &o) const {return frac(p*o.q-q*o.p,q*o.q);}
    frac operator *(const frac &o) const {return frac(p*o.p,q*o.q);}

    bool operator ==(const frac &o) const {return p*o.q==q*o.p;}
    bool operator !=(const frac &o) const {return p*o.q==q*o.p;}
    bool operator <(const frac &o) const {return p*o.q<q*o.p;}
    bool operator >(const frac &o) const {return p*o.q>q*o.p;}
    bool operator <=(const frac &o) const {return p*o.q<=q*o.p;}
    bool operator >=(const frac &o) const {return p*o.q>=q*o.p;}

    double get_val() const {return 1.0*p/q;}

    friend istream& operator >> (istream& is, frac &b){
        is >> b.p >> b.q;
        return is;
    }
    friend ostream& operator << (ostream& os, const frac &b){
        return os << b.p << '/' << b.q;
    }
};
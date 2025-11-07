#pragma once

#include "library/mod/modint.hpp"
#include "library/poly/convolution.hpp"

template<typename mint>
struct poly: vector<mint>{
    using vector<mint>::vector;

    poly(const vector<mint> &v):vector<mint>(v){}

    poly operator += (const poly &o){
        if(o.size()>this->size()) this->resize(o.size());
        for(int i=0; i<(int)o.size(); ++i) (*this)[i]+=o[i];
        return *this;
    }
    poly operator -= (const poly &o){
        if(o.size()>this->size()) this->resize(o.size());
        for(int i=0; i<(int)o.size(); ++i) (*this)[i]-=o[i];
        return *this;
    }
    poly operator *= (const poly &o){
        return *this=convolution<mint>(*this,o);
    }
    poly operator + (const poly &o) const {return poly(*this)+=o;}
    poly operator - (const poly &o) const {return poly(*this)-=o;}
    poly operator * (const poly &o) const {return poly(*this)*=o;}
};
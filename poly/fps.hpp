#pragma once

#include "library/mod/modint_basic.hpp"
#include "library/poly/convolution.hpp"

template<typename mint>
struct fps: vector<mint>{
    using vector<mint>::vector;

    fps(const vector<mint> &v):vector<mint>(v){}

    fps operator += (const fps &o){
        if(o.size()>this->size()) this->resize(o.size());
        for(int i=0; i<(int)o.size(); ++i) (*this)[i]+=o[i];
        return *this;
    }
    fps operator += (const mint &o){
        if(this->empty()) this->resize(1);
        (*this)[0]+=o;
        return *this;
    }
    fps operator -= (const fps &o){
        if(o.size()>this->size()) this->resize(o.size());
        for(int i=0; i<(int)o.size(); ++i) (*this)[i]-=o[i];
        return *this;
    }
    fps operator -= (const mint &o){
        if(this->empty()) this->resize(1);
        (*this)[0]-=o;
        return *this;
    }
    fps operator *= (const fps &o){
        return *this=convolution<mint>(*this,o);
    }
    fps operator *= (const mint &o){
        for(int i=0; i<(int)this->size(); ++i) (*this)[i]*=o;
        return *this;
    }
    fps operator + (const fps &o) const {return fps(*this)+=o;}
    fps operator + (const mint &o) const {return fps(*this)+=o;}
    fps operator - (const fps &o) const {return fps(*this)-=o;}
    fps operator - (const mint &o) const {return fps(*this)-=o;}
    fps operator * (const fps &o) const {return fps(*this)*=o;}
    fps operator * (const mint &o) const {return fps(*this)*=o;}

    fps operator - () const {return fps()-*this;}

    fps topos(){
        while((this->size())&&(this->back())==mint(0)){
            this->pop_back();
        }
        return *this;
    }

    friend ostream& operator << (ostream& os, const fps &P){
        int n=P.size();
        for(int i=0; i<n; ++i){
            os << P[i];
            if(i+1<n) os << ' ';
        }
        return os;
    }
};
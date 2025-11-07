#pragma once

#include "library/mod/modint.hpp"
#include "library/mod/modint_basic.hpp"
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
    poly operator += (const mint &o){
        if(this->empty()) this->resize(1);
        (*this)[0]+=o;
        return *this;
    }
    poly operator -= (const poly &o){
        if(o.size()>this->size()) this->resize(o.size());
        for(int i=0; i<(int)o.size(); ++i) (*this)[i]-=o[i];
        return *this;
    }
    poly operator -= (const mint &o){
        if(this->empty()) this->resize(1);
        (*this)[0]-=o;
        return *this;
    }
    poly operator *= (const poly &o){
        return *this=convolution<mint>(*this,o);
    }
    poly operator *= (const mint &o){
        for(int i=0; i<(int)this->size(); ++i) (*this)[i]*=o;
        return *this;
    }
    poly operator + (const poly &o) const {return poly(*this)+=o;}
    poly operator + (const mint &o) const {return poly(*this)+=o;}
    poly operator - (const poly &o) const {return poly(*this)-=o;}
    poly operator - (const mint &o) const {return poly(*this)-=o;}
    poly operator * (const poly &o) const {return poly(*this)*=o;}
    poly operator * (const mint &o) const {return poly(*this)*=o;}

    poly interval(int l, int r){
        assert(l<=r&&r<=(int)this->size());
        poly res(this->begin()+l,this->begin()+r);
        return res;
    }

    poly inverse(){
        int n=this->size();
        assert((*this)[0]!=0);
        poly res(1,(*this)[0].inv());
        poly b=*this;
        for(int m=1; m<n; m<<=1){
            if(n<m*2) b.resize(m*2);
            poly v1=b.interval(0,m*2),v2=res;
            v1*=v2;
            v1.resize(m*2);
            v1*=v2;
            res.resize(m*2);
            for(int i=0; i<m; ++i) res[i]+=res[i];
            for(int i=0; i<m*2; ++i) res[i]-=v1[i];
        }
        res.resize(n);
        return res;
    }

    poly derivative(){
        int n=this->size();
        poly res(n-1);
        for(int i=0; i<n-1; ++i) res[i]=(*this)[i+1]*mint(i+1);
        return res;
    }

    poly integral(){
        int n=this->size();
        poly res(n+1);
        for(int i=0; i<n; ++i) res[i+1]=(*this)[i]*(inv<mint>(i+1));
        return res;
    }
};
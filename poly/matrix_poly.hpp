#pragma once

#include "library/poly/fps.hpp"
#include "library/poly/fps_basic.hpp"

template<typename fps>
struct matrix_poly: vector<fps>{
    using vector<fps>::vector;
    using mint=typename fps::value_type;
    static_assert(mint::can_ntt());

    inline static NTT<mint> ntt;

    int h,w;

    matrix_poly(){
        h=w=0;
    }

    void set_hw(int _h, int _w){
        h=_h,w=_w;
        this->resize(h*w);
        for(int i=0; i<h*w; ++i) (*this)[i]=fps();
    }
    static matrix_poly id(int n){
        matrix_poly res;
        res.set_hw(n,n);
        for(int i=0; i<n*n; ++i) res[i]=fps();
        for(int i=0; i<n; ++i) res[i*n+i]=fps({1});
        return res;
    }
    matrix_poly transpose(){
        matrix_poly res;
        res.set_hw(w,h);
        for(int i=0; i<h; ++i) for(int j=0; j<w; ++j){
            res[j*h+i]=(*this)[i*w+j];
        }
        return res;
    }
    matrix_poly operator += (const matrix_poly &mat2){
        assert(h==mat2.h&&w==mat2.w);
        for(int i=0; i<h*w; ++i) (*this)[i]+=mat2[i];
        return *this;
    }
    matrix_poly operator -= (const matrix_poly &mat2){
        assert(h==mat2.h&&w==mat2.w);
        for(int i=0; i<h*w; ++i) (*this)[i]-=mat2[i];
        return *this;
    }
    matrix_poly operator *= (matrix_poly mat2){
        assert(w==mat2.h);
        mat2=mat2.transpose();
        matrix_poly res;
        res.set_hw(h,mat2.h);
        matrix_poly mat1=*this;
        int mx_n=0,mx_m=0;
        for(int i=0; i<h*w; ++i){
            mx_n=max(mx_n,(int)mat1[i].size());
        }
        for(int i=0; i<mat2.h*w; ++i){
            mx_m=max(mx_m,(int)mat2[i].size());
        }
        int t=0;
        while((1<<t)<mx_n+mx_m-1) t++;
        for(int i=0; i<h*w; ++i){
            mat1[i].resize(1<<t);
        }
        for(int i=0; i<mat2.h*w; ++i){
            mat2[i].resize(1<<t);
        }
        for(int i=0; i<h*mat2.h; ++i){
            res[i].resize(1<<t);
        }
        for(int i=0; i<h*w; ++i){
            ntt.dft(mat1[i],t);
        }
        for(int i=0; i<mat2.h*w; ++i){
            ntt.dft(mat2[i],t);
        }
        for(int i=0; i<h; ++i) for(int j=0; j<mat2.h; ++j){
            for(int k=0; k<w; ++k){
                res[i*mat2.h+j]+=dot<mint>(mat1[i*w+k],mat2[j*w+k]);
            }
        }
        for(int i=0; i<h*mat2.h; ++i){
            ntt.dft(res[i],t,true);
            res[i].topos();
        }
        *this=res;
        return *this;
    }
    matrix_poly operator *= (const mint &x){
        for(int i=0; i<h*w; ++i) (*this)[i]*=x;
        return *this;
    }
    matrix_poly operator + (const matrix_poly &o) const {return matrix_poly(*this)+=o;}
    matrix_poly operator - (const matrix_poly &o) const {return matrix_poly(*this)-=o;}
    matrix_poly operator * (const matrix_poly &o) const {return matrix_poly(*this)*=o;}
    matrix_poly operator * (const mint &o) const {return matrix_poly(*this)*=o;}
    matrix_poly operator ^(ll x){
        assert(h==w);
        assert(x>=0);
        matrix_poly res=matrix_poly::id(h),cur=*this;
        for(; x; cur*=cur,x>>=1) if(x&1) res*=cur;
        return res;
    }
};
#pragma once

#include "library/simd/matmul.hpp"

template<typename T>
struct matrix: vector<T>{
    using vector<T>::vector;

    inline static bool is_modint=false;
    int h,w;

    matrix(){
        h=w=0;
    }

    void set_hw(int _h, int _w){
        h=_h,w=_w;
        this->resize(h*w);
        for(int i=0; i<h*w; ++i) (*this)[i]=0;
    }
    static matrix id(int n){
        matrix res;
        res.set_hw(n,n);
        for(int i=0; i<n*n; ++i) res[i]=0;
        for(int i=0; i<n; ++i) res[i*n+i]=1;
        return res;
    }
    matrix transpose(){
        matrix res;
        res.set_hw(w,h);
        for(int i=0; i<h; ++i) for(int j=0; j<w; ++j){
            res[j*h+i]=(*this)[i*w+j];
        }
        return res;
    }
    matrix operator += (const matrix &mat2){
        assert(h==mat2.h&&w==mat2.w);
        for(int i=0; i<h*w; ++i) (*this)[i]+=mat2[i];
        return *this;
    }
    matrix operator -= (const matrix &mat2){
        assert(h==mat2.h&&w==mat2.w);
        for(int i=0; i<h*w; ++i) (*this)[i]-=mat2[i];
        return *this;
    }
    matrix operator *= (matrix mat2){
        assert(w==mat2.h);
        mat2=mat2.transpose();
        matrix res;
        res.set_hw(h,mat2.h);
        bool simd=0;
        #ifndef i_am_noob
        if(is_modint) simd=1;
        #endif
        if(simd){
            vector<int> vec1(h*w),vec2(mat2.h*w);
            for(int i=0; i<h*w; ++i){
                vec1[i]=(*this)[i].x;
            }
            for(int i=0; i<mat2.h*w; ++i){
                vec2[i]=mat2[i].x;
            }
            vector<int> vec3=matmul_mod(h,mat2.h,w,vec1,vec2,T::get_mod());
            for(int i=0; i<h*mat2.h; ++i){
                res[i]=vec3[i];
            }
        }
        else{
            for(int i=0; i<h; ++i) for(int j=0; j<mat2.h; ++j){
                for(int k=0; k<w; ++k){
                    res[i*mat2.h+j]+=(*this)[i*w+k]*mat2[j*w+k];
                }
            }
        }
        *this=res;
        return *this;
    }
    matrix operator *= (const T &x){
        for(int i=0; i<h*w; ++i) (*this)[i]*=x;
        return *this;
    }
    matrix operator + (const matrix &o) const {return matrix(*this)+=o;}
    matrix operator - (const matrix &o) const {return matrix(*this)-=o;}
    matrix operator * (const matrix &o) const {return matrix(*this)*=o;}
    matrix operator * (const T &o) const {return matrix(*this)*=o;}
    matrix operator ^(ll x){
        assert(h==w);
        assert(x>=0);
        matrix res=matrix::id(h),cur=*this;
        for(; x; cur*=cur,x>>=1) if(x&1) res*=cur;
        return res;
    }
    friend ostream& operator << (ostream& os, const matrix &x){
        for(int i=0; i<x.h; ++i){
            for(int j=0; j<x.w; ++j){
                os << x[i*x.w+j];
                if(j+1==x.w) os << '\n';
                else os << ' ';
            }
        }
        return os;
    }
};
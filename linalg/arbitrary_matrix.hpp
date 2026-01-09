#pragma once

template<typename T>
struct matrix: vector<T>{
    using vector<T>::vector;

    int h,w;

    matrix(){
        h=w=0;
    }

    void set_hw(int _h, int _w){
        h=_h,w=_w;
        this->resize(h*w);
        for(int i=0; i<h*w; ++i) (*this)[i]=0;
    }
    /*
    static matrix id(int n){
        matrix res;
        res.set_hw(n,n);
        for(int i=0; i<n*n; ++i) res[i]=0;
        for(int i=0; i<n; ++i) res[i*n+i]=1;
        return res;
    }
    */
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
    matrix operator *= (const matrix &mat2){
        assert(w==mat2.h);
        matrix res;
        res.set_hw(h,mat2.w);
        for(int i=0; i<h; ++i) for(int j=0; j<mat2.w; ++j){
            for(int k=0; k<w; ++k) res[i*mat2.w+j]+=(*this)[i*w+k]*mat2[k*mat2.w+j];
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
};
#pragma once

template<typename T, int N>
struct matrix: array<array<T,N>,N>{
    matrix(){for(int i=0; i<N; ++i) for(int j=0; j<N; ++j) (*this)[i][j]=0;}
    void init(vector<T> vec){
        for(int i=0; i<N; ++i) for(int j=0; j<N; ++j) (*this)[i][j]=vec[i*N+j];
    }
    static constexpr matrix id(){
        matrix res;
        for(int i=0; i<N; ++i) res[i][i]=1;
        return res;
    }
    matrix t(){
        matrix res;
        for(int i=0; i<N; ++i) for(int j=0; j<N; ++j) res[i][j]=(*this)[j][i];
        return res;
    }
    matrix operator += (const matrix &mat2){
        for(int i=0; i<N; ++i) for(int j=0; j<N; ++j) (*this)[i][j]+=mat2[i][j];
        return *this;
    }
    matrix operator -= (const matrix &mat2){
        for(int i=0; i<N; ++i) for(int j=0; j<N; ++j) (*this)[i][j]-=mat2[i][j];
        return *this;
    }
    matrix operator *= (matrix mat2){
        mat2=mat2.t();
        matrix res;
        for(int i=0; i<N; ++i) for(int j=0; j<N; ++j){
            for(int k=0; k<N; ++k) res[i][j]+=(*this)[i][k]*mat2[j][k];
        }
        *this=res;
        return *this;
    }
    matrix operator + (const matrix &o) const {return matrix(*this)+=o;}
    matrix operator - (const matrix &o) const {return matrix(*this)-=o;}
    matrix operator * (const matrix &o) const {return matrix(*this)*=o;}
    matrix operator ^(ll x){
        assert(x>=0);
        matrix res=id(),cur=*this;
        for(; x; cur*=cur,x>>=1) if(x&1) res*=cur;
        return res;
    }
    friend ostream& operator << (ostream& os, const matrix &x){
        for(int i=0; i<N; ++i){
            for(int j=0; j<N; ++j){
                os << x[i][j];
                if(j+1==N) os << '\n';
                else os << ' ';
            }
        }
        return os;
    }
};

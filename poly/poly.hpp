#pragma once

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

    poly operator - () const {return poly()-*this;}

    poly interval(int l, int r){
        assert(l<=r&&r<=(int)this->size());
        poly res(this->begin()+l,this->begin()+r);
        return res;
    }

    poly topos(){
        while((this->size())&&(this->back())==mint(0)){
            this->pop_back();
        }
        return *this;
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

    pair<poly,poly> divide(poly b){
        // return {quotient, remainder}
        poly a=*this;
        int n=a.size(),m=b.size(),k=n-m+1;
        if(n<m) return {poly({0}),a};
        if(mint::ntt_data().first<0||m<=50){
            poly q(k),r;
            mint tmp=b[m-1].inv();
            for(int i=k-1; i>=0; --i){
                q[i]=a[m-1+i]*tmp;
                for(int j=0; j<m; ++j){
                    a[i+j]-=q[i]*b[j];
                }
                assert(a[m-1+i]==mint(0));
            }
            a.resize(m-1);
            r=a;
            return {q,r};
        }
        poly ra=a,rb=b;
        reverse(all(ra)),reverse(all(rb));
        ra.resize(k),rb.resize(k);
        poly q=ra*rb.inverse();
        q.resize(k);
        reverse(all(q));
        poly r=a-b*q;
        r.resize(m-1);
        return {q,r};
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

    poly ln(){
        // a[0] = 1
        int n=this->size();
        if(n==1) return poly({0});
        poly d=derivative();
        poly b=*this;
        b.pop_back();
        poly res=d*b.inverse();
        res.resize(n-1);
        return res.integral();
    }

    poly exp(){
        // a[0] = 0
        int n=this->size();
        poly q(1,1);
        poly b=*this;
        b[0]+=1;
        for(int m=1; m<n; m<<=1){
            if(n<m*2) b.resize(m*2);
            poly g=b.interval(0,m*2),h=q;
            h.resize(m*2),h=h.ln();
            g-=h;
            q*=g;
            q.resize(m*2);
        }
        q.resize(n);
        return q;
    }

    poly pow_naive(ll k){
        int n=this->size();
        poly b=*this,res={1};
        for(; k; b*=b,k>>=1,b.resize(n)) if(k&1) res*=b,res.resize(n);
        return res;
    }

    int low(){
        int n=this->size(),m=0;
        while(m<n&&(*this)[m]==0) m++;
        if(m>=n) return -1;
        return m;
    }
    poly shift(int n){
        poly res(n,0);
        res.insert(res.end(),this->begin(),this->end());
        return res;
    }

    poly pow(ll k){ // 0^0 = 1
        int n=this->size();
        if(k==0){
            poly res(n);
            return res[0]=1,res;
        }
        int m=low();
        if(m){
            if(m==-1||k>=n||k*m>=n) return poly(n);
            int lft=n-k*m;
            poly b=interval(m,m+lft);
            b=b.pow(k);
            b=b.shift(k*m);
            return b;
        }
        poly b=*this;
        mint base=b[0].pow(k),inv=b[0].inv();
        b*=inv;
        b=b.ln();
        if(b.empty()) b.pb(0);
        b*=k;
        b=b.exp();
        b*=base;
        return b;
    }

    poly pow_sparse(int k, int n){ // 0^0 = 1
        if(k==0){
            poly res(n);
            return res[0]=1,res;
        }
        int t=this->size(),m=low();
        if(m){
            if(m==-1||k>=n||1ll*k*m>=n) return poly(n);
            int lft=n-k*m;
            poly b=interval(m,t);
            b=b.pow_sparse(k,lft);
            b=b.shift(k*m);
            return b;
        }
        poly res(n,0);
        res[0]=(*this)[0].pow(k);
        mint inv_a0=(*this)[0].inv();
        for(int i=1; i<n; ++i){
            for(int j=1; j<t; ++j){
                if(i-j>=0) res[i]-=res[i-j]*(i-j)*(*this)[j];
            }
            for(int j=1; j<t; ++j){
                if(i-j>=0) res[i]+=res[i-j]*(*this)[j]*j*k;
            }
            res[i]*=inv_a0*inv<mint>(i);
        }
        return res;
    }

    friend ostream& operator << (ostream& os, const poly &P){
        int n=P.size();
        for(int i=0; i<n; ++i){
            os << P[i];
            if(i+1<n) os << ' ';
        }
        return os;
    }
};
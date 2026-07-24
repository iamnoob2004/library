#pragma once

#include "nt/basic.hpp"
#include "nt/array_of_floor.hpp"

template<typename T>
struct dirichlet{
    ll n;
    int m;
    array_of_floor arr;
    vector<T> val;

    dirichlet(ll _n=0){init(_n);}

    void init(ll _n){
        n=_n;
        arr.init(n);
        m=arr.m;
        val.assign(m,T());
    }

    dirichlet operator *= (const dirichlet &o){
        assert(n==o.n);
        vector<T> f=val,g=o.val,h(m,T());
        auto op=[&](int xl, int xr, int yl, int yr, int zl, int zr){
            T tmp=(f[xr]-f[xl-1])*(g[yr]-g[yl-1]);
            h[zl-1]+=tmp,h[zr]-=tmp;
        };
        for(int i=1; i<m; ++i){
            // x biggest (y,z < x), z = m-i
            for(int j=1; j<m; ++j){
                int ii=m-arr.get_id(arr[j]*arr[m-i]);
                if(ii<=max(j,m-i)) break;
                op(max(j,m-i)+1,ii,j,j,m-i,m-i);
            }
            // z biggest (x,y <= z), x = i
            // y biggest (x <= y, z < y), x = i
            for(int j=1; j<m; ++j){
                int k=m-arr.get_id(arr[i]*arr[j]);
                if(k<max(i,j)) break;
                op(i,i,j,j,max(i,j),k);
                op(i,i,max(i,j+1),k,j,j);
            }
        }
        for(int i=m-1; i>=1; --i) h[i]=h[i-1];
        h[0]=0;
        for(int i=1; i<m; ++i) h[i]+=h[i-1];
        reverse(h.begin()+1,h.end());
        val=h;
        return *this;
    }

    dirichlet operator /= (const dirichlet &o){
        assert(n==o.n);
        vector<T> f(m,T()),g=o.val,h=val;
        reverse(h.begin()+1,h.end());
        for(int i=0; i<m-1; ++i){
            h[i]=h[i+1]-h[i];
        }
        h[m-1]=-h[m-1];
        auto op=[&](int xl, int xr, int yl, int yr, int zl, int zr){
            T tmp=(f[xr]-f[xl-1])*(g[yr]-g[yl-1]);
            h[zl-1]-=tmp,h[zr]+=tmp;
        };
        T c=g[1],cinv=T(1)/c;
        for(int i=1; i<m; ++i){
            for(int j=2; j<m; ++j){
                int ii=m-arr.get_id(arr[j]*arr[m-i]);
                if(ii<=max(j,m-i)) break;
                op(max(j,m-i)+1,ii,j,j,m-i,m-i);
            }
            // determine f[i]
            if(arr[i]<=n/arr[i]){
                // (f[i]-f[i-1])*g[1]+h[m-i] = 0
                f[i]=-h[m-i]*cinv+f[i-1];
            }
            else{
                f[i]=-h[m-i]*cinv+f[m-i];
            }
            if(n/(arr[1]*arr[m-i])>arr[m-i]){
                op(m-i+1,m-arr.get_id(arr[1]*arr[m-i]),1,1,m-i,m-i);
            }
            for(int j=1; j<m; ++j){
                int k=m-arr.get_id(arr[i]*arr[j]);
                if(k<max(i,j)) break;
                op(i,i,j,j,max(i,j),k);
                op(i,i,max(i,j+1),k,j,j);
            }
            assert(h[m-i]==mint(0));
        }
        val=f;
        return *this;
    }

    dirichlet operator * (const dirichlet &o) const {return dirichlet(*this)*=o;}
    dirichlet operator / (const dirichlet &o) const {return dirichlet(*this)/=o;}
};

template<typename T>
dirichlet<T> dirichlet_identity(ll n){
    dirichlet<T> res(n);
    for(int i=1; i<res.m; ++i) res.val[i]=T(1);
    return res;
}

template<typename T>
dirichlet<T> dirichlet_one(ll n){
    dirichlet<T> res(n);
    for(int i=1; i<res.m; ++i) res.val[i]=T(res.arr[i]);
    return res;
}

template<typename T>
dirichlet<T> dirichlet_n(ll n){
    dirichlet<T> res(n);
    for(int i=1; i<res.m; ++i) res.val[i]=sum_of_power_small<T>(res.arr[i],2);
    return res;
}

template<typename T>
dirichlet<T> dirichlet_phi(ll n){
    return dirichlet_n<T>(n)/dirichlet_one<T>(n);
}

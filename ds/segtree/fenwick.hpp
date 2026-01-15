#pragma once

#include "library/algebra/monoid/add.hpp"

// abelian group
template<typename group>
struct fenwick{
    using T=typename group::value_type;

    int n;
    vector<T> vec,tree;
    T tot;

    fenwick(){build(0);}
    fenwick(int _n){build(_n);}
    fenwick(int _n, T val){
        build(n,val);
    }
    fenwick(vector<T> _vec){
        build(_vec);
    }

    void build(int _n, T val=group::id()){
        n=_n;
        vector<T> _vec;
        _vec.reserve(n);
        _vec.assign(n,val);
        build(_vec);
    }
    void build(vector<T> _vec){
        n=_vec.size();
        vec.reserve(n);
        vec.assign(n,group::id());
        tree.reserve(n+1);
        tree.assign(n+1,group::id());
        tot=group::id();
        for(int i=0; i<n; ++i){
            add(i,_vec[i]);
        }
    }

    void set(int p, T x){
        assert(0<=p&&p<n);
        add(p,group::op(x,group::inv(vec[p])));
    }
    void add(int p, T x){
        assert(0<=p&&p<n);
        vec[p]=group::op(vec[p],x);
        tot=group::op(tot,x);
        for(++p; p<=n; p+=p&-p) tree[p]=group::op(tree[p],x);
    }

    T sum(int r){
        assert(0<=r&&r<=n);
        T res=group::id();
        for(; r; r-=r&-r) res=group::op(res,tree[r]);
        return res;
    }
    T sum(int l, int r){
        assert(l<=r);
        return group::op(sum(r),group::inv(sum(l)));
    }
};
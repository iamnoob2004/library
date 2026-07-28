#pragma once

template<typename monoid_action>
struct lazy_segtree{
    using value_type=monoid_action;

    using monoid_S=typename monoid_action::monoid_S;
    using monoid_F=typename monoid_action::monoid_F;
    using S=typename monoid_action::S;
    using F=typename monoid_action::F;

    lazy_segtree():lazy_segtree(0){}
    lazy_segtree(int _n):lazy_segtree(vector<S>(_n,monoid_S::id())){}

    int n,m,log;

    vector<S> vec,tree;
    vector<F> tag;

    lazy_segtree(vector<S> _vec){
        n=_vec.size();
        log=0;
        while(n>(1<<log)) log++;
        m=1<<log;
        vec=_vec;

        tree.assign(m*2,monoid_S::id());
        tag.assign(m*2,monoid_F::id());
        for(int i=0; i<n; ++i){
            tree[m+i]=vec[i];
        }

        for(int i=m-1; i>=1; --i){
            tree[i]=monoid_S::op(tree[i*2],tree[i*2+1]);
        }
    }

    void pull(int k){tree[k]=monoid_S::op(tree[k*2],tree[k*2+1]);}
    void give_tag(int k, F f){
        tree[k]=monoid_action::act(f,tree[k]);
        tag[k]=monoid_F::op(f,tag[k]);
    }
    void push(int k){
        give_tag(k*2,tag[k]);
        give_tag(k*2+1,tag[k]);
        tag[k]=monoid_F::id();
    }

    void set(int p, S x){
        p+=m;
        for(int i=log; i>=1; --i) push(p>>i);
        tree[p]=x;
        for(int i=1; i<=log; ++i) pull(p>>i);
    }
    S get(int p){
        p+=m;
        for(int i=log; i>=1; --i) push(p>>i);
        return tree[p];
    }

    S prod(int l, int r){
        if(l==r) return monoid_S::id();
        l+=m,r+=m;
        
        for(int i=log; i>=1; --i){
            if(((l>>i)<<i)!=l) push(l>>i);
            if(((r>>i)<<i)!=r) push(r>>i);
        }

        S resl=monoid_S::id(),resr=monoid_S::id();

        while(l<r){
            if(l&1) resl=monoid_S::op(resl,tree[l++]);
            if(r&1) resr=monoid_S::op(tree[--r],resr);
            l>>=1,r>>=1;
        }

        return monoid_S::op(resl,resr);
    }

    S all_prod(){return tree[1];}

    void apply(int p, F f){
        p+=m;
        for(int i=log; i>=1; --i) push(p>>i);
        tree[p]=monoid_action::act(f,tree[p]);
        for(int i=1; i<=log; ++i) pull(p>>i);
    }

    void apply(int l, int r, F f){
        if(l==r) return;
        l+=m,r+=m;
        
        for(int i=log; i>=1; --i){
            if(((l>>i)<<i)!=l) push(l>>i);
            if(((r>>i)<<i)!=r) push(r>>i);
        }

        {
            int l2=l,r2=r;
            while(l<r){
                if(l&1) give_tag(l++,f);
                if(r&1) give_tag(--r,f);
                l>>=1,r>>=1;
            }
            l=l2,r=r2;
        }

        for(int i=1; i<=log; ++i){
            if(((l>>i)<<i)!=l) pull(l>>i);
            if(((r>>i)<<i)!=r) pull(r>>i);
        }
    }
};
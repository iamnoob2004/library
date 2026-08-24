#pragma once

template<typename monoid, bool persistent=false>
struct dynamic_segtree{
    using value_type=monoid;

    using S=typename monoid::S;
    using F=function<S(ll,ll)>;
    using node_id=int;

    #define MID ((l+r)>>1)

    struct node{
        node_id ch[2];
        S x;
    };

    F default_prod;
    vector<node> pool;
    ll l_rt,r_rt;
    node_id rt=-1;

    node_id new_node(S x){
        pool.pb({{-1,-1},x});
        return (int)pool.size()-1;
    }
    node_id copy_node(node_id k){
        if(!persistent) return k;
        pool.pb(pool[k]);
        return (int)pool.size()-1;
    }

    dynamic_segtree():dynamic_segtree(0,0){}
    dynamic_segtree(ll l, ll r, F _default_prod=[](ll, ll){return monoid::id();}){
        l_rt=l,r_rt=r;
        default_prod=_default_prod;
    }
    dynamic_segtree(vector<S> _vec):dynamic_segtree(0,_vec.size()){
        if(!_vec.empty()) rt=build_rec(l_rt,r_rt,_vec);
    }
    node_id build_rec(ll l, ll r, vector<S>& vec){
        if(l+1==r) return new_node(vec[l]);
        ll mid=MID;
        node_id lc=build_rec(l,mid,vec),rc=build_rec(mid,r,vec);
        node_id res=new_node(monoid::op(pool[lc].x,pool[rc].x));
        pool[res].ch[0]=lc,pool[res].ch[1]=rc;
        return res;
    }

    S get_val(node_id k, ll l, ll r){
        return k==-1?default_prod(l,r):pool[k].x;
    }

    void pull(int k, ll l, ll r){
        ll mid=MID;
        pool[k].x=monoid::op(get_val(pool[k].ch[0],l,mid),get_val(pool[k].ch[1],mid,r));
    }

    // returns new root
    node_id set(node_id k, ll p, S x){
        assert(l_rt<=p&&p<r_rt);
        return set_rec(k,l_rt,r_rt,p,x);
    }
    node_id set_rec(node_id k, ll l, ll r, ll p, S x){
        if(k==-1) k=new_node(default_prod(l,r));
        else k=copy_node(k);
        if(l+1==r){
            pool[k].x=x;
            return k;
        }
        ll mid=MID;
        int branch=p>=mid;
        node_id ch=set_rec(pool[k].ch[branch],branch?mid:l,branch?r:mid,p,x);
        pool[k].ch[branch]=ch;
        pull(k,l,r);
        return k;
    }

    S prod(node_id k, ll ql, ll qr){
        assert(l_rt<=ql&&ql<=qr&&qr<=r_rt);
        return prod_rec(k,l_rt,r_rt,ql,qr);
    }
    S prod_rec(node_id k, ll l, ll r, ll ql, ll qr){
        if(ql>=r||qr<=l) return monoid::id();
        if(k==-1) return default_prod(max(l,ql),min(r,qr));
        if(ql<=l&&qr>=r) return pool[k].x;
        ll mid=MID;
        return monoid::op(
            prod_rec(pool[k].ch[0],l,mid,ql,qr),
            prod_rec(pool[k].ch[1],mid,r,ql,qr));
    }

    #undef MID
};
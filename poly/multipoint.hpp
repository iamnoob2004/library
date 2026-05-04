#pragma once

#include "library/poly/fps_inv.hpp"
#include "library/poly/middle_product.hpp"

template<typename mint>
struct subproduct_tree{
    int m,k;
    vector<vector<mint>> tree;
    subproduct_tree(vector<mint> vec){
        m=vec.size(),k=1;
        while(k<m) k<<=1;
        tree.resize(k*2);
        for(int i=0; i<k; ++i){
            // reversed
            if(i<m) tree[k+i]={1,-vec[i]};
            else tree[k+i]={1,0};
        }
        for(int i=k-1; i>=1; --i){
            tree[i]=convolution<mint>(tree[i<<1],tree[i<<1|1]);
        }
    }

    // transposition principle
    vector<mint> eval(vector<mint> f){
        int n=f.size();
        if(n==0) return vector<mint>(m,mint(0));
        f.resize(n*2-1);
        vector<vector<mint>> tree2(k*2);
        tree2[1]=tree[1];
        tree2[1].resize(n);
        tree2[1]=fps_inv<mint>(tree2[1]);
        tree2[1]=middle_product<mint>(f,tree2[1]);
        tree2[1].resize(k);
        for(int i=1; i<k; ++i){
            tree2[i<<1]=middle_product<mint>(tree2[i],tree[i<<1|1]);
            tree2[i<<1|1]=middle_product<mint>(tree2[i],tree[i<<1]);
        }
        vector<mint> res(m);
        for(int i=0; i<m; ++i) res[i]=tree2[k+i][0];
        return res;
    }

    vector<mint> interpolation(vector<mint> val){
        assert(m==(int)val.size());
        vector<mint> a(m);
        for(int i=0; i<m; ++i) a[i]=tree[1][m-1-i]*(i+1);
        a=eval(a);

        vector<vector<mint>> tree2(k*2);
        for(int i=0; i<k; ++i) tree2[k+i]={(i<m?val[i]/a[i]:0)};
        for(int i=k-1; i>=1; --i){
            tree2[i]=convolution<mint>(tree2[i<<1],tree[i<<1|1]);
            auto tmp=convolution<mint>(tree2[i<<1|1],tree[i<<1]);
            for(int j=0; j<(int)tree2[i].size(); ++j) tree2[i][j]+=tmp[j];
        }
        tree2[1].resize(m);
        reverse(tree2[1].begin(),tree2[1].end());
        return tree2[1];
    }
};

template<typename mint>
vector<mint> multipoint_eval(vector<mint> f, vector<mint> x){
    if(x.empty()) return {};
    subproduct_tree<mint> tree(x);
    return tree.eval(f);
}

template<typename mint>
vector<mint> multipoint_interpolation(vector<mint> x, vector<mint> y){
    if(x.empty()) return {};
    subproduct_tree<mint> tree(x);
    return tree.interpolation(y);
}

// calculate f(ar^k) for 0 <= k < m
template<typename mint>
vector<mint> multipoint_eval_on_geom_seq(vector<mint> f, mint a, mint r, int m){
    int n=f.size();
    if(m==0) return {};
    auto eval=[&](mint x) -> mint{
        mint res,pw=1;
        for(int i=0; i<n; ++i){
            res+=f[i]*pw;
            pw*=x;
        }
        return res;
    };
    if(r==mint(0)){
    	vector<mint> res(m);
    	for(int i=1; i<m; ++i) res[i]=f[0];
    	res[0]=eval(a);
    	return res;
    }
    if(n<49||m<49){
    	vector<mint> res(m);
    	for(int i=0; i<m; ++i,a*=r) res[i]=eval(a);
    	return res;
    }
    {
    	mint pw=1;
    	for(int i=0; i<n; ++i,pw*=a) f[i]*=pw;
    }
    auto calc=[&](mint q, int k){
        // calc q^binom(i,2) for 0 <= i < k
        vector<mint> res(k);
        mint pw=1;
        res[0]=1;
        for(int i=0; i<k-1; ++i,pw*=q){
            res[i+1]=res[i]*pw;
        }
        return res;
    };
    vector<mint> vec1=calc(r.inv(),max(n,m)),vec2=calc(r,n+m);
    vector<mint> c=f;
    for(int i=0; i<n; ++i){
        c[i]=f[i]*vec1[i];
    }
    c=middle_product<mint>(vec2,c);
    vector<mint> res(m);
    for(int i=0; i<m; ++i){
        res[i]=c[i]*vec1[i];
    }
    return res;
}
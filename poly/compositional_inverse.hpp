#pragma once

#include "library/poly/poly.hpp"
#include "library/poly/power_projection.hpp"

template<typename mint>
poly<mint> compositional_inverse(poly<mint> f){
    int n=((int)f.size())-1;
    if(n==-1) return {};
    assert(f[0]==mint(0));
    if(n==0) return f;
    assert(f[1]!=mint(0));
    mint c=f[1],inv_c=c.inv();
    for(int i=0; i<=n; ++i) f[i]*=inv_c;
    /*
    f(g(x)) = x
    Lagrange inversion formula:
    n[x^n]f(x)^i = i[x^-i]g(x)^(-n) = i[x^(n-i)](g(x)/x)^(-n)
    if [x^1]g = 1, we can calculate g(x)/x from (g(x)/x)^(-n) by poly pow
    it suffices to calculate n[x^n]f(x)^i for i=1,2,...,n
    */
    vector<mint> w(n+1);
    w[n]=1;
    vector<mint> val=power_projection<mint>(f,w,n);
    poly<mint> g(n);
    for(int i=1; i<=n; ++i) g[n-i]=val[i]*mint(n)*inv<mint>(i);
    g=g.pow((-inv<mint>(n)).x);
    g.insert(g.begin(),mint(0));
    mint pow_inv_c=1;
    for(int i=0; i<=n; ++i) g[i]*=pow_inv_c,pow_inv_c*=inv_c;
    return g;
}
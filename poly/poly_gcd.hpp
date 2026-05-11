#pragma once

#include "library/poly/matrix_poly.hpp"
#include "library/poly/poly_div.hpp"

namespace half_gcd{

template<typename mint>
using matrix=matrix_poly<fps<mint>>;

template<typename mint>
matrix<mint> step(vector<mint> q){
    matrix<mint> res;
    res.set_hw(2,2);
    res[1]=res[2]={1};
    res[3]-=q;
    return res;
}

// go until one poly is smaller then deg(a[0])/2
template<typename mint>
matrix<mint> hgcd(matrix<mint> a){
    assert(a.h==2&&a.w==1);
    assert(a[0].size()>a[1].size()&&a[1].size()>0);
    int m=a[0].size()/2;
    if((int)a[1].size()<=m) return matrix<mint>::id(2);
    matrix<mint> a_right;
    a_right.set_hw(2,1);
    for(int i=0; i<2; ++i){
        a_right[i]=vector<mint>(a[i].begin()+m,a[i].end());
    }
    auto M1=hgcd(a_right);
    a=M1*a;
    if((int)a[1].size()<=m) return M1;
    matrix<mint> Q=step<mint>(poly_div<mint>(a[0],a[1]).first);
    M1=Q*M1,a=Q*a;
    if((int)a[1].size()<=m) return M1;
    int m2=2*m+1-(int)a[0].size();
    for(int i=0; i<2; ++i){
        a_right[i]=vector<mint>(a[i].begin()+m2,a[i].end());
    }
    auto M2=hgcd(a_right);
    return M2*M1;
}

template<typename mint>
matrix<mint> full_gcd(matrix<mint> a){
    assert(a.h==2&&a.w==1);
    assert(a[0].size()>a[1].size()&&a[1].size()>0);
    auto M1=hgcd<mint>(a);
    a=M1*a;
    if(a[1].empty()) return M1;
    matrix<mint> Q=step<mint>(poly_div<mint>(a[0],a[1]).first);
    M1=Q*M1,a=Q*a;
    if(a[1].empty()) return M1;
    return full_gcd(a)*M1;
}

// ax + by = gcd(a,b), {gcd(a,b),x,y}
template<typename mint>
array<vector<mint>,3> poly_extgcd(vector<mint> a, vector<mint> b){
    while(!a.empty()&&a.back()==mint(0)) a.pop_back();
    while(!b.empty()&&b.back()==mint(0)) b.pop_back();
    if(a.empty()) return {b,{},{1}};
    if(b.empty()) return {a,{1},{}};
    matrix<mint> Q=step<mint>(poly_div<mint>(a,b).first);
    auto M=Q;
    matrix<mint> vec;
    vec.set_hw(2,1);
    vec[0]=a,vec[1]=b;
    vec=Q*vec;
    if(!vec[1].empty()) M=full_gcd<mint>(vec)*M;
    fps<mint> _a(a),_b(b);
    fps<mint> g=M[0]*_a+M[1]*_b;
    g.topos();
    return {g,M[0],M[1]};
}

}// namespace half_gcd

using half_gcd::poly_extgcd;
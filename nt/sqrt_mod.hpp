#pragma once

template<typename mint>
mint legendre(mint x){
    int p=mint::get_mod();
    return x.pow((p-1)/2);
}

template<typename mint>
bool is_quadratic_residue(mint x){
    int p=mint::get_mod();
    if(p==2||x==mint(0)) return true;
    return legendre<mint>(x)==mint(1);
}

template<typename mint>
mint sqrt_mod(mint x){
    int p=mint::get_mod();
    if(p==2||x==mint(0)) return x;
    assert(is_quadratic_residue(x));
    mt19937 rng(49);
    mint t;
    while(is_quadratic_residue(t*t-x)) t=mint(rng()%p);
    // bostan-mori
    mint f0(1),f1(0),g1=t*mint(-2),g2=x;
    for(int i=(p+1)/2; i; i>>=1){
        mint nf0,nf1,ng1,ng2;
        if(i&1){
            nf0=f1-f0*g1,nf1=f1*g2;
        }
        else{
            nf0=f0,nf1=f0*g2-f1*g1;
        }
        ng1=g2*mint(2)-g1*g1,ng2=g2*g2;
        f0=nf0,f1=nf1,g1=ng1,g2=ng2;
    }
    return f0;
}
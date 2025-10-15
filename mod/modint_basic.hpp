#pragma once

template<typename mint>
mint fac(int n){
    static const int mod=mint::get_mod();
    static vector<mint> res={1,1};
    if(n>=mod) return 0;
    assert(n>=0);
    while(n>=(int)res.size()) res.push_back(res.back()*res.size());
    return res[n];
}

template<typename mint>
mint inv(int n){
    static const int mod=mint::get_mod();
    static vector<mint> res={0,1};
    assert(n>=0&&n<mod);
    while(n>=(int)res.size()) res.push_back(res[mod%res.size()]*(mod-mod/res.size()));
    return res[n];
}

template<typename mint>
mint ifac(int n){
    static const int mod=mint::get_mod();
    static vector<mint> res={1,1};
    if(n>=mod) return 0;
    assert(n>=0);
    while(n>=(int)res.size()) res.push_back(res.back()*inv<mint>(res.size()));
    return res[n];
}

template<typename mint>
mint C(int n, int m){
    if(m<0||m>n) return 0;
    return fac<mint>(n)*ifac<mint>(m)*ifac<mint>(n-m);
}
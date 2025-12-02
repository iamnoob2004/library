#pragma once

// ax + by = gcd(a,b), {gcd(a,b),x,y}
template<typename T>
array<T,3> extgcd(T a, T b){
    T x1=1,y1=0,x2=0,y2=1;
    while(b!=0){
        T q=a/b;
        a%=b;
        swap(a,b);
        T x3=x1-x2*q,y3=y1-y2*q;
        x1=x2,y1=y2,x2=x3,y2=y3;
    }
    return {a,x1,y1};
}

template<typename T>
T modinv(T x, T m){
    auto [g,val1,val2]=extgcd<T>(x,m);
    assert(g==1);
    if(val1<0) val1+=m;
    return val1;
}
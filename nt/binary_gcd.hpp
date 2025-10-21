#pragma once

template<typename T>
T bgcd(T a, T b){
    if(a==0) return b;
    if(b==0) return a;
    int az=__builtin_ctzll(a);
    int bz=__builtin_ctzll(b);
    int shift=min(az,bz);
    b>>=bz;
    while(a!=0){
        a>>=az;
        T diff=b-a;
        az=__builtin_ctzll(diff);
        b=min(a,b);
        a=abs(diff);
    }
    return b<<shift;
}
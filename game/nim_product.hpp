#pragma once

struct nim_product_solver_64{
    using u16=uint16_t;
    using u64=uint64_t;

    static const int small=1<<16;
    static const u16 g=-1;
    u16 Log[small],pw[small];

    template<bool flag=1>
    u64 nim_mul(u64 a, u64 b, int mid=32){
        if(a<=1||b<=1) return a*b;
        if(flag&&a<small&&b<small){
            int x=Log[a]+Log[b];
            if(x>=small-1) x-=small-1;
            return pw[x];
        }
        u64 s=(1ULL<<mid)-1;
        u64 a0=a&s,a1=a>>mid,b0=b&s,b1=b>>mid;
        u64 val1=nim_mul<flag>(a0,b0,mid/2),val2=nim_mul<flag>(a0^a1,b0^b1,mid/2),val3=nim_mul<flag>(a1,b1,mid/2);
        return val1^((val1^val2)<<mid)^nim_mul<flag>(val3,1ULL<<(mid-1),mid/2);
    }

    nim_product_solver_64(){
        pw[0]=1;
        for(int i=1; i<small-1; ++i) pw[i]=nim_mul<0>(pw[i-1],g,8);
        for(int i=0; i<small-1; ++i) Log[pw[i]]=i;
    }

    u64 get(u64 a, u64 b){return nim_mul(a,b);}
};


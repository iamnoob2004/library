#pragma once

template<typename mint>
struct NTT{
    static constexpr int m=mint::get_mod(),N=mint::ntt_data().first,g=mint::ntt_data().second;

    mint w[N+1];

    NTT(){
        w[N]=g;
        for(int i=N-1; i>=0; --i) w[i]=w[i+1]*w[i+1];
    }

    void trans(vector<mint> &a, int k, bool inv=false){
        for(int i=1,j=0; i<(1<<k); ++i){
            for(int t=1<<(k-1); (j^=t)<t; t>>=1);
            if(i<j) swap(a[i],a[j]);
        }
        for(int L=1,step=2; L<=k; ++L,step<<=1){
            for(int i=0; i<(1<<k); i+=step){
                mint cur(1),dw=w[L];
                for(int j=i,j2=i+(step>>1); j<i+(step>>1); ++j,++j2,cur*=dw){
                    mint tmp=a[j2]*cur;
                    a[j2]=a[j]-tmp;
                    a[j]+=tmp;
                }
            }
        }
        if(inv){
            reverse(a.begin()+1,a.end());
            mint inv=mint(1<<k).inv();
            for(int i=0; i<(1<<k); ++i) a[i]*=inv;
        }
    }
};
#pragma once

template<typename mint>
struct NTT{
    static constexpr int mod=mint::get_mod(),N=mint::ntt_data().first,g=mint::ntt_data().second;

    mint w[2][N+1],w1[2][N],w2[2][N];
    using u64=uint64_t;

    NTT(){
        w[0][N]=g,w[1][N]=w[0][N].inv();
        for(int i=N-1; i>=0; --i){
            w[0][i]=w[0][i+1]*w[0][i+1];
            w[1][i]=w[1][i+1]*w[1][i+1];
        }
        mint prod=1,iprod=1;
        for(int i=0; i<N-1; ++i){
            w1[0][i]=w[0][i+2]*prod;
            w1[1][i]=w[1][i+2]*iprod;
            prod*=w[1][i+2];
            iprod*=w[0][i+2];
        }
        prod=iprod=1;
        for(int i=0; i<N-2; ++i){
            w2[0][i]=w[0][i+3]*prod;
            w2[1][i]=w[1][i+3]*iprod;
            prod*=w[1][i+3];
            iprod*=w[0][i+3];
        }
    }

    void dft(vector<mint> &a, int k, bool inv=false){
        assert((int)a.size()==(1<<k));
        int n=1<<k;
        if(!inv){
            int len=0;
            while(len<k){
                if(k-len==1){
                    int p=1<<(k-len-1);
                    mint rot=1;
                    for(int s=0; s<(1<<len); ++s){
                        int offset=s<<(k-len);
                        for(int i=0; i<p; ++i){
                            mint l=a[i+offset],r=a[i+offset+p]*rot;
                            a[i+offset]+=r,a[i+offset+p]=l-r;
                        }
                        rot*=w1[0][__lg(~s&-~s)];
                    }
                    len++;
                }
                else{
                    int p=1<<(k-len-2);
                    mint rot=1,imag=w[0][2];
                    u64 mod2=((u64)mod)*mod;
                    for(int s=0; s<(1<<len); ++s){
                        mint rot2=rot*rot,rot3=rot2*rot;
                        int offset=s<<(k-len);
                        for(int i=0; i<p; ++i){
                            u64 a0=a[i+offset].x,a1=u64(a[i+offset+p].x)*rot.x,a2=u64(a[i+offset+2*p].x)*rot2.x,a3=u64(a[i+offset+3*p].x)*rot3.x;
                            u64 tmp=(a1+mod2-a3)%mod*imag.x;
                            u64 na2=mod2-a2;
                            a[i+offset]=a0+a2+a1+a3;
                            a[i+offset+p]=a0+a2+(2*mod2-(a1+a3));
                            a[i+offset+2*p]=a0+na2+tmp;
                            a[i+offset+3*p]=a0+na2+(mod2-tmp);
                        }
                        rot*=w2[0][__lg(~s&-~s)];
                    }
                    len+=2;
                }
            }
        }
        else{
            int len=k;
            while(len){
                if(len==1){
                    len--;
                    int p=1<<(k-len-1);
                    mint rot=1;
                    for(int s=0; s<(1<<len); ++s){
                        int offset=s<<(k-len);
                        for(int i=0; i<p; ++i){
                            u64 l=a[i+offset].x,r=a[i+offset+p].x;
                            a[i+offset]=l+r,a[i+offset+p]=(mod+l-r)*rot.x;
                        }
                        rot*=w1[1][__lg(~s&-~s)];
                    }
                }
                else{
                    len-=2;
                    int p=1<<(k-len-2);
                    mint rot=1,imag=w[1][2];
                    for(int s=0; s<(1<<len); ++s){
                        mint rot2=rot*rot,rot3=rot2*rot;
                        int offset=s<<(k-len);
                        for(int i=0; i<p; ++i){
                            u64 a0=a[i+offset].x,a1=a[i+offset+p].x,a2=a[i+offset+2*p].x,a3=a[i+offset+3*p].x;
                            u64 tmp=(mod+a2-a3)*imag.x%mod;
                            a[i+offset]=a0+a1+a2+a3;
                            a[i+offset+p]=(a0+mod-a1+tmp)*rot.x;
                            a[i+offset+2*p]=(a0+a1+(2*mod-(a2+a3)))*rot2.x;
                            a[i+offset+3*p]=(a0+2*mod-a1-tmp)*rot3.x;
                        }
                        rot*=w2[1][__lg(~s&-~s)];
                    }
                }
            }
            mint inv=mint(n).inv();
            for(int i=0; i<n; ++i) a[i]*=inv;
        }
    }

    void transpose_dft(vector<mint> &a, int k, bool inv=false){
        int n=1<<k;
        if(!inv){
            dft(a,k,1);
            reverse(a.begin()+1,a.end());
            for(auto &x: a) x*=n;
        }
        else{
            reverse(a.begin()+1,a.end());
            dft(a,k,0);
            mint inv_n=mint(n).inv();
            for(auto &x: a) x*=inv_n;
        }
    }
};
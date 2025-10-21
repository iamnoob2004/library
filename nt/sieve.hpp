#pragma once

template<int N>
struct sieve{
    int lpd[N],pr[N],n_primes;
    sieve(){
        memset(lpd,-1,sizeof lpd);
        n_primes=0;
        for(int i=2; i<N; ++i){
            if(lpd[i]==-1) lpd[i]=i,pr[n_primes++]=i;
            for(int j=0; j<n_primes; ++j){
                int &p=pr[j];
                if(p*i>=N) break;
                lpd[p*i]=p;
                if(i%p==0) break;
            }
        }
    }
    bool isp(int n){return n==lpd[n];}
    vector<int> factorize(int n){
        vector<int> res;
        while(n>1){
            res.pb(lpd[n]);
            n/=lpd[n];
        }
        return res;
    }
    int divisor_sum(int n){
        vector<int> vec=factorize(n);
        int res=1,cur=1,q;
        for(int i=0; i<vec.size(); ++i){
            if(i==0||vec[i]!=vec[i-1]) res*=cur,cur=vec[i]+1,q=vec[i];
            else q*=vec[i],cur+=q;
        }
        res*=cur;
        return res;
    }
    int totient(int n){
        vector<int> vec=factorize(n);
        int res=n;
        for(int i=0; i<vec.size(); ++i){
            if(i==0||vec[i]!=vec[i-1]) res=res/vec[i]*(vec[i]-1);
        }
        return res;
    }
};
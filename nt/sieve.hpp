#pragma once

struct sieve{
    int N;
    vector<int> lpd,pr;
    void init(int _N){
        N=_N;
        lpd.assign(N,-1);
        for(int i=2; i<N; ++i){
            if(lpd[i]==-1) lpd[i]=i,pr.pb(i);
            for(int p: pr){
                if(p*i>=N) break;
                lpd[p*i]=p;
                if(i%p==0) break;
            }
        }
    }
    sieve(){}
    sieve(int _N){
        init(_N);
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
        for(int i=0; i<(int)vec.size(); ++i){
            if(i==0||vec[i]!=vec[i-1]) res*=cur,cur=vec[i]+1,q=vec[i];
            else q*=vec[i],cur+=q;
        }
        res*=cur;
        return res;
    }
    int totient(int n){
        vector<int> vec=factorize(n);
        int res=n;
        for(int i=0; i<(int)vec.size(); ++i){
            if(i==0||vec[i]!=vec[i-1]) res=res/vec[i]*(vec[i]-1);
        }
        return res;
    }
};
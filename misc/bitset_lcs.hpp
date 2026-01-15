#pragma once

#include "library/ds/dynamic_bitset.hpp"

int bitset_lcs(vector<int> a, vector<int> b, int mx){
    mx++;
    int n=a.size();
    vector<dynamic_bitset> vec(mx,dynamic_bitset(n+1));
    for(int i=0; i<n; ++i) vec[a[i]].set(i);

    dynamic_bitset res(n+1),tmp(n+1);
    for(int y: b){
        tmp=vec[y];
        tmp|=res;
        res<<=1;
        res.set(0);
        res=tmp-res;
        res.flip();
        res&=tmp;
    }
    return res.count();
}

vector<int> bitset_lcs_recover(vector<int> a, vector<int> b, int mx){
    mx++;
    int n=a.size(),m=b.size();
    int k=sqrt(m);
    k=max(k,1);
    vector<dynamic_bitset> vec(mx,dynamic_bitset(n+1));
    for(int i=0; i<n; ++i) vec[a[i]].set(i);

    vector<dynamic_bitset> res;
    dynamic_bitset cur(n+1),tmp(n+1);
    for(int j=0; j<=m; ++j){
        if(j%k==0){
            res.pb(cur);
        }
        if(j==m) break;
        tmp=vec[b[j]];
        tmp|=cur;
        cur<<=1;
        cur.set(0);
        cur=tmp-cur;
        cur.flip();
        cur&=tmp;
    }
    int ans_len=cur.count();
    vector<int> lcs;
    int cur_pos=n;
    for(int t=((int)res.size())-1; t>=0; --t){
        cur=res[t];
        vector<dynamic_bitset> vec2;
        for(int j=t*k; j<(t+1)*k&&j<m; ++j){
            tmp=vec[b[j]];
            tmp|=cur;
            cur<<=1;
            cur.set(0);
            cur=tmp-cur;
            cur.flip();
            cur&=tmp;
            vec2.pb(cur);
        }
        for(int j=((int)vec2.size())-1; j>=0; --j){
            cur_pos=vec2[j].find_prev(cur_pos);
            if(cur_pos>=0&&a[cur_pos]==b[t*k+j]){
                lcs.pb(a[cur_pos]);
                cur_pos--;
            }
        }
    }
    reverse(lcs.begin(),lcs.end());
    assert(ans_len==lcs.size());
    return lcs;
}
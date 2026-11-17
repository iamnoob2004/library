#pragma once

struct dsu_undo{
    vector<int> lead,siz;
    dsu_undo(int n=0){
        lead.resize(n);
        siz.resize(n);
        for(int i=0; i<n; ++i) lead[i]=i,siz[i]=1;
    }
    int leader(int x){
        while(lead[x]!=x) x=lead[x];
        return x;
    }
    bool same(int x, int y){return leader(x)==leader(y);}
    int size(int x){return siz[leader(x)];}
    struct data{
        int x,par,siz;
        data(int _x, int _par, int _siz):x(_x),par(_par),siz(_siz){}
    };
    vector<data> stk;
    int merge(int x, int y){
        x=leader(x),y=leader(y);
        if(x==y) return x;
        if(siz[x]<siz[y]) swap(x,y);
        stk.pb(data(y,x,siz[y]));
        lead[y]=x;
        siz[x]+=siz[y];
        return x;
    }
    int undo(){
        if(stk.empty()) return -1;
        siz[stk.back().par]-=stk.back().siz;
        lead[stk.back().x]=stk.back().x;
        int res=stk.back().x;
        stk.pop_back();
        return res;
    }
};
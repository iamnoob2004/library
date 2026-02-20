#pragma once

#include "library/graph/graph.hpp"

template<typename graph_type>
struct Tree{
    int n,maxlog;
    graph_type g;

    vector<vector<int>> anc;
    vector<int> par,tl,tr,dep,ord,siz,ch,head;

    Tree(int _n=0){
        n=g.n=_n;
    }
    void read(bool weighted=false, int offset=1){
        g.read_tree(weighted,offset);
    }

    bool gen_anc; // maxlog, anc
    bool gen_hld; // ch, head
    int rt;

    // ch and head don't need much space, so default true
    void build(int _rt=0, bool need_hld=1, bool need_anc=0){
        g.build();
        rt=_rt;
        gen_hld=need_hld;
        gen_anc=need_anc;
        maxlog=0;
        while((1<<maxlog)<=n) maxlog++;
        par.resize(n);
        tl.resize(n),tr.resize(n);
        dep.resize(n);
        siz.resize(n);
        ord.clear();
        if(need_anc){
            anc.resize(n);
            for(int i=0; i<n; ++i){
                anc[i].resize(maxlog);
            }
        }
        if(need_hld){
            ch.resize(n);
            head.resize(n);
            for(int i=0; i<n; ++i){
                ch[i]=-1;
            }
        }
        _time=-1;
        dfs1(rt,-1);
        dfs2(rt,-1);
        if(need_hld){
            head[rt]=rt;
            for(int u: ord){
                for(auto &e: g.get_adj(u)) if(e.to!=par[u]){
                    if(e.to==ch[u]) head[e.to]=head[u];
                    else head[e.to]=e.to;
                }
            }
        }
    }

    void dfs1(int u, int fa){
        par[u]=fa;
        siz[u]=1;
        if(fa==-1) dep[u]=0;
        else dep[u]=dep[fa]+1;
        for(auto &e: g.get_adj(u)) if(e.to!=fa){
            dfs1(e.to,u);
            siz[u]+=siz[e.to];
            if(gen_hld){
                if(ch[u]==-1||siz[e.to]>siz[ch[u]]) ch[u]=e.to;
            }
        }
    }

    int _time;
    void dfs2(int u, int fa){
        ord.push_back(u);
        tl[u]=++_time;
        if(gen_anc){
            anc[u][0]=fa;
            for(int i=1; i<maxlog; ++i){
                if(anc[u][i-1]==-1) anc[u][i]=-1;
                else anc[u][i]=anc[anc[u][i-1]][i-1];
            }
        }
        if(gen_hld&&ch[u]!=-1){
            dfs2(ch[u],u);
        }
        for(auto &e: g.get_adj(u)) if(e.to!=fa&&(!gen_hld||e.to!=ch[u])){
            dfs2(e.to,u);
        }
        tr[u]=_time;
    }

    bool is_anc(int u, int v){return tl[u]<=tl[v]&&tr[u]>=tr[v];}
    int get_anc(int u, int x){
        assert(gen_anc);
        for(int i=maxlog-1; i>=0; --i) if(u!=-1&&(x>>i&1)) u=anc[u][i];
        return u;
    }
    int lca(int u, int v){
        assert(gen_anc);
        if(is_anc(u,v)) return u;
        for(int i=maxlog-1; i>=0; --i) if(anc[u][i]!=-1&&!is_anc(anc[u][i],v)) u=anc[u][i];
        return par[u];
    }

    // returns {vertex list, edge list}
    // each edge stores {parent, child}
    pair<vector<int>,vector<pii>> virtual_tree(vector<int> vec){
        sort(all(vec),[&](int u, int v){return tl[u]<tl[v];});
        int s=vec.size();
        for(int i=0; i<s-1; ++i) vec.push_back(lca(vec[i],vec[i+1]));
        sort(all(vec),[&](int u, int v){return tl[u]<tl[v];});
        vec.resize(unique(all(vec))-vec.begin());
        vector<pii> edges; // {parent, child}
        vector<int> stk;
        for(int u: vec){
            while(!stk.empty()&&!is_anc(stk.back(),u)) stk.pop_back();
            if(!stk.empty()) edges.pb({stk.back(),u});
            stk.push_back(u);
        }
        return {vec,edges};
    }
};
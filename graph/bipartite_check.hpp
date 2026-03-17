#pragma once

#include "library/graph/graph.hpp"

template<typename graph_type>
pair<vector<int>,vector<int>> bipartite_check(graph_type &g){
    int n=g.n;
    assert(g.gen);
    vector<int> vis(n),col(n);
    bool ok=1;
    for(int i=0; i<n; ++i) if(!vis[i]){
        if(!ok) break;
        queue<int> q;
        q.push(i),vis[i]=1;
        while(!q.empty()&&ok){
            int u=q.front(); q.pop();
            for(auto &e: g.get_adj(u)){
                if(!ok) break;
                if(vis[e.to]){
                    if(col[u]==col[e.to]) ok=0;
                    continue;
                }
                col[e.to]=col[u]^1;
                q.push(e.to),vis[e.to]=1;
            }
        }
    }
    if(!ok) return {{-1},{-1}};
    vector<int> res[2];
    for(int i=0; i<n; ++i){
    	res[col[i]].pb(i);
    }
    return {res[0],res[1]};
}
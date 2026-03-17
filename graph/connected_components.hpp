#pragma once

#include "library/graph/graph.hpp"

// undirected only
// return component id for each vertex
template<typename graph_type>
vector<int> connected_components(graph_type &g){
    static_assert(!graph_type::directed_type);
    int n=g.n;
    assert(g.gen);
    vector<int> vis(n),res(n);
    int cur_id=0;
    for(int i=0; i<n; ++i) if(!vis[i]){
        queue<int> q;
        q.push(i),vis[i]=1;
        while(!q.empty()){
            int u=q.front(); q.pop();
            res[u]=cur_id;
            for(auto &e: g.get_adj(u)) if(!vis[e.to]){
                q.push(e.to),vis[e.to]=1;
            }
        }
        cur_id++;
    }
    return res;
}

template<typename graph_type>
vector<graph_type> connected_components_subgraph(graph_type &g){
    int n=g.n;
    if(n==0) return {};
    vector<int> component_id=connected_components<graph_type>(g);
    int k=*max_element(component_id.begin(),component_id.end());
    k++;
    vector<int> component_size(k),id2(n);
    for(int i=0; i<n; ++i){
        id2[i]=component_size[component_id[i]]++;
    }
    vector<graph_type> res(k);
    for(int i=0; i<k; ++i){
        res[i]=graph(component_size[i]);
    }
    for(auto &e: g.edges){
        int i=component_id[e.from];
        res[i].add_edge(id2[e.from],id2[e.to],e.cost,e.id);
    }
    for(int i=0; i<k; ++i){
        res[i].build();
    }
    return res;
}
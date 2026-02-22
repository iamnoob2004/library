#pragma once

template<typename T>
struct edge{
    int from,to;
    T cost;
    int id;
};

template<typename T=int,bool directed=false>
struct graph{
    using weight_type=T;
    int n,m;
    vector<edge<T>> edges,edge_list;
    vector<int> start_ptr;
    bool gen;

    void add_edge(int u, int v, T w=T(1), int id=-1){
        assert(0<=u&&u<n&&0<=v&&v<n&&!gen);
        if(id==-1) id=m;
        edges.pb(edge<T>{u,v,w,id});
        m++;
    }
    graph():n(0),m(0),gen(0){}
    graph(int _n):n(_n),m(0),gen(0){}
    void read_graph(int _m, bool weighted=false, int offset=1){
        assert(!gen);
        for(int i=0; i<_m; ++i){
            int u,v; cin >> u >> v;
            u-=offset,v-=offset;
            if(!weighted) add_edge(u,v);
            else{
                T w; cin >> w;
                add_edge(u,v,w);
            }
        }
    }
    void read_tree(bool weighted=false, int offset=1){
        read_graph(n-1,weighted,offset);
    }

    void build(){
        assert(!gen);
        gen=1;
        start_ptr.resize(n+1);
        for(auto &e: edges){
            start_ptr[e.from+1]++;
            if(!directed) start_ptr[e.to+1]++;
        }
        for(int i=0; i<n; ++i) start_ptr[i+1]+=start_ptr[i];
        auto cur=start_ptr;
        edge_list.resize(cur.back());
        for(auto e: edges){
            edge_list[cur[e.from]++]=e;
            if(!directed){
                swap(e.from,e.to);
                edge_list[cur[e.from]++]=e;
            }
        }
    }

    vector<edge<T>> get_adj(int u){
        assert(gen);
        vector<edge<T>> vec(edge_list.begin()+start_ptr[u],edge_list.begin()+start_ptr[u+1]);
        return vec;
    }
};
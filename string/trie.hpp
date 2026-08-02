#pragma once

template<typename T=int, int C=26>
struct trie{
    vector<array<int,C>> ch;
    vector<T> vec;
    int nw_node(){
        int res=ch.size();
        array<int,C> arr;
        for(int i=0; i<C; ++i) arr[i]=-1;
        ch.pb(arr);
        vec.pb(T());
        return res;
    }
    void init(){
        ch.clear();
        vec.clear();
        nw_node();
    }
    trie(){
        init();
    }
    int insert(vector<int> a){
        int cur=0;
        for(int i: a){
            if(ch[cur][i]==-1) ch[cur][i]=nw_node();
            cur=ch[cur][i];
        }
        return cur;
    }
    int insert(string s, char first){
        vector<int> a;
        for(auto c: s) a.pb(c-first);
        return insert(a);
    }

    vector<int> dfs_order(){
        ord.clear();
        internal_dfs(0);
        return ord;
    }
  private:
    vector<int> ord;
    void internal_dfs(int x){
        ord.pb(x);
        for(int i=0; i<C; ++i) if(ch[x][i]!=-1){
            internal_dfs(ch[x][i]);
        }
    }
};
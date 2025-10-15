#pragma once

struct grid{
    static constexpr int dx[8]={-1,0,1,0,-1,-1,1,1},dy[8]={0,-1,0,1,-1,1,-1,1};
    int n,m;

    void set_size(int _n, int _m){
        n=_n,m=_m;
    }

    grid(){}
    grid(int _n, int _m){
        set_size(_n,_m);
    }

    bool inside(int i, int j){return i>=0&&i<n&&j>=0&&j<m;}
    int id(int i, int j){return i*m+j;}
    pair<int,int> pos(int x){return {x/m,x%m};}

    bool is_adj4(int u, int v){
        if(u<0||u>=n*m||v<0||v>=n*m) return false;
        auto [x1,y1]=pos(u);
        auto [x2,y2]=pos(v);
        return abs(x1-x2)+abs(y1-y2)==1;
    }

    vector<int> adj(int u, int d=4){
        auto [x,y]=pos(u);
        vector<int> res;
        for(int i=0; i<d; ++i){
            int nx=x+dx[i],ny=y+dy[i];
            if(inside(nx,ny)) res.pb(id(nx,ny));
        }
        return res;
    }

    int area() const {
        return n*m;
    }
};
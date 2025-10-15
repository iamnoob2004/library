#pragma once

struct dsu{
	int n;
	vector<int> par,siz;

	void build(int _n){
		n=_n;
		par.resize(n);
		siz.resize(n);
		for(int i=0; i<n; ++i){
			par[i]=i;
			siz[i]=1;
		}
	}

	dsu(){}
	dsu(int _n){
		build(_n);
	}

	int find(int x){
		if(par[x]==x) return x;
		return par[x]=find(par[x]);
	}
	int size(int x){
		return siz[find(x)];
	}
	bool merge(int x, int y){
		x=find(x),y=find(y);
		if(x==y) return false;
		if(siz[x]<siz[y]) swap(x,y);
		siz[x]+=siz[y];
		par[y]=x;
		return true;
	}
};
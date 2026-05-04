#pragma once

template<typename fps>
fps all_product(vector<fps> vec){
	if(vec.empty()) return fps({1});
    auto dfs=[&](auto &self, int l, int r) -> fps{
        if(l+1==r) return vec[l];
        int mid=(l+r)>>1;
        return self(self,l,mid)*self(self,mid,r);
    };
    return dfs(dfs,0,vec.size());
}
#pragma once

template<typename poly>
poly all_product(vector<poly> vec){
    auto dfs=[&](auto &self, int l, int r) -> poly{
        if(l+1==r) return vec[l];
        int mid=(l+r)>>1;
        return self(self,l,mid)*self(self,mid,r);
    };
    return dfs(dfs,0,vec.size());
}
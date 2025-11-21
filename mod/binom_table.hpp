#pragma once

template<typename T>
vector<vector<T>> binom_table(int n){
	vector<vector<T>> res(n+1);
	for(int i=0; i<=n; ++i){
		res[i].resize(i+1);
		for(int j=0; j<=i; ++j){
			if(j==0||j==i) res[i][j]=T(1);
			else res[i][j]=res[i-1][j-1]+res[i-1][j];
		}
	}
	return res;
}
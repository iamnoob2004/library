#pragma once

template<typename monoid>
struct sliding_window_aggregation{
    using T=typename monoid::value_type;
    using value_type=T;

    int siz;
    vector<T> rval,lprod;
    T rprod;

    sliding_window_aggregation(){
        siz=0;
        rval={};
        lprod={monoid::id()};
        rprod=monoid::id();
    }

    void push(T x){
        rprod=monoid::op(rprod,x);
        rval.pb(x);
        siz++;
    }

    void pop(){
        assert(siz>0);
        lprod.pop_back();
        if(lprod.empty()){
            lprod={monoid::id()};
            while(!rval.empty()){
                lprod.pb(monoid::op(rval.back(),lprod.back()));
                rval.pop_back();
            }
            rprod=monoid::id();
            lprod.pop_back();
        }
        siz--;
    }

    T prod(){
        return monoid::op(lprod.back(),rprod);
    }
};
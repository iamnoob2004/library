#pragma once

#include "nt/bigint.hpp"

struct bigint_binary {
    using Int=bigint_binary;
    using u32=uint32_t;
    using u64=uint64_t;
    static constexpr u32 base=1<<30;
    static constexpr int base_digits=30;

    vector<u32> a;
    bool neg=false;

    bigint_binary(){}
    bigint_binary(ll x){
        a.clear();
        neg=x<0;
        u64 y=neg?u64(-(x+1))+1:u64(x);
        while(y){
            a.push_back(y%base);
            y/=base;
        }
    }

    Int& shrink(){
        while(!a.empty()&&a.back()==0) a.pop_back();
        if(a.empty()) neg=false;
        return *this;
    }

    friend istream& operator >> (istream& is, Int &b){
        bigint_1e9 tmp;
        is >> tmp;
        b=Int();
        for(int i=(int)tmp.a.size()-1; i>=0; --i){
            b*=Int(bigint_1e9::base);
            b+=Int(tmp.a[i]);
        }
        b.neg=tmp.neg;
        return is;
    }

    friend ostream& operator << (ostream& os, const Int &b){
        bigint_1e9 tmp;
        for(int i=(int)b.a.size()-1; i>=0; --i){
            tmp*=bigint_1e9(Int::base);
            tmp+=bigint_1e9(b.a[i]);
        }
        tmp.neg=b.neg;
        return os << tmp;
    }

    explicit operator bool() const {
        return !a.empty();
    }
    explicit operator int() const {
        u64 x=0;
        for(int i=0; i<(int)a.size()&&i<2; ++i){
            x=x*base+a[i];
        }
        if(neg&&x==1ll<<31){
            return (-1)<<31;
        }
        x&=(1ll<<31)-1;
        if(neg) x=-x;
        return x;
    }

    // a < b: -1, a = b: 0, a > b: 1
    friend int compare(const Int &a, const Int &b){
        if(a.neg!=b.neg) return a.neg?-1:1;
        return compare_abs(a,b)*(a.neg?-1:1);
    }

    friend int operator ==(const Int &a, const Int &b){
        return a.neg==b.neg&&a.a==b.a;
    }
    friend int operator <(const Int &a, const Int &b){
        return compare(a,b)<0;
    }
    friend int operator >(const Int &a, const Int &b){return b<a;}
    friend int operator !=(const Int &a, const Int &b){return !(a==b);}
    friend int operator <=(const Int &a, const Int &b){return !(b<a);}
    friend int operator >=(const Int &a, const Int &b){return !(a<b);}

    Int& operator &= (const Int &o){
        if(a.size()>o.a.size()) a.resize(o.a.size());
        for(int i=0; i<(int)a.size(); ++i) a[i]&=o.a[i];
        return shrink();
    }
    friend Int& operator & (Int a, const Int &b){return a&=b;}

    Int& operator >>= (const int &o){
        assert(o>=0);
        vector<u32> res(a.size());
        for(int i=0; i<(int)a.size(); ++i){
            for(int j=0; j<base_digits; ++j){
                if(a[i]>>j&1){
                    int nw=i*base_digits+j-o;
                    if(nw>=0){
                        res[nw/base_digits]|=1u<<(nw%base_digits);
                    }
                }
            }
        }
        a=res;
        return shrink();
    }

    Int& operator - (){
        if(a.size()) neg=!neg;
        return *this;
    }

    Int& operator += (const Int& o){
        if(neg==o.neg){
            add_abs(o);
            return *this;
        }
        if(compare_abs(*this,o)>=0){
            sub_abs(o);
        }
        else{
            Int tmp=o;
            tmp.sub_abs(*this);
            *this=tmp;
        }
        return shrink();
    }
    Int& operator -= (const Int &o){
        Int tmp=o;
        return *this+=(-tmp);
    }

    Int& operator *= (const Int &o){
        if(a.empty()||o.a.empty()) return *this=bigint_binary();
        vector<u32> res(a.size()+o.a.size());
        for(int i=0; i<(int)a.size(); ++i){
            u64 x=0,t=a[i];
            for(int j=0; j<(int)o.a.size()||x; ++j){
                u64 nw=res[i+j]+x;
                if(j<(int)o.a.size()) nw+=t*o.a[j];
                res[i+j]=nw%base;
                x=nw/base;
            }
        }
        a=res;
        neg^=o.neg;
        return shrink();
    }

private:
    // ignores sign
    void add_abs(const Int& o){
        u32 x=0;
        if(a.size()<o.a.size()){
            a.resize(o.a.size());
        }
        for(int i=0; i<(int)a.size(); ++i){
            if(i<(int)o.a.size()) a[i]+=o.a[i];
            a[i]+=x;
            x=0;
            if(a[i]>=base) a[i]-=base,x++;
        }
        if(x) a.push_back(x);
        shrink();
    }
    // |*this| >= |o|
    void sub_abs(const Int& o){
        u32 x=0;
        for(int i=0; i<(int)a.size(); ++i){
            int y=a[i];
            if(i<(int)o.a.size()) y-=o.a[i];
            y-=x;
            x=0;
            if(y<0) y+=base,x++;
            a[i]=y;
        }
        shrink();
    }
    // a < b: -1, a = b: 0, a > b: 1
    friend int compare_abs(const Int &a, const Int &b){
        if(a.a.size()!=b.a.size()) return a.a.size()<b.a.size()?-1:1;
        for(int i=((int)a.a.size())-1; i>=0; --i){
            if(a.a[i]!=b.a[i]){
                return a.a[i]<b.a[i]?-1:1;
            }
        }
        return 0;
    }
};

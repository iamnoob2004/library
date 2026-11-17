#pragma once

struct bigint_1e9 {
    using Int=bigint_1e9;
    using u32=uint32_t;
    using u64=uint64_t;
    static constexpr u32 base=1000000000;
    static constexpr int base_digits=9;

    vector<u32> a;
    bool neg=false;

    bigint_1e9(){}
    bigint_1e9(ll x){
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
        string s;
        is >> s;
        b=Int();
        reverse(s.begin(),s.end());
        if(s.back()=='-') b.neg=true,s.pop_back();
        int j=0,x=1;
        for(int i=0; i<(int)s.size(); ++i){
            assert(s[i]>='0'&&s[i]<='9');
            while(j>=(int)b.a.size()) b.a.push_back(0);
            b.a[j]+=(s[i]-'0')*x;
            x*=10;
            if(x==base){
                j++,x=1;
            }
        }
        return is;
    }

    friend ostream& operator << (ostream& os, const Int &b){
        string s;
        if(b.a.empty()){
            s="0";
            return os << s;
        }
        for(int i=0; i<(int)b.a.size(); ++i){
            u32 x=b.a[i];
            for(int j=0; j<base_digits; ++j){
                s+='0'+(x%10);
                x/=10;
            }
        }
        while(!s.empty()&&s.back()=='0') s.pop_back();
        if(b.neg) s+='-';
        reverse(s.begin(),s.end());
        return os << s;
    }
    
    explicit operator bool() const {
        return !a.empty();
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
        if(a.empty()||o.a.empty()) return *this=bigint_1e9();
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

#pragma once

template<typename T=double>
struct surreal{
    static constexpr T inf=T(1ll<<61);

    T x;
    surreal():x(0){}
    surreal(T _x):x(_x){}

    bool operator == (const surreal &o) const {return x==o.x;}
    bool operator != (const surreal &o) const {return x!=o.x;}
    bool operator < (const surreal &o) const {return x<o.x;}
    bool operator <= (const surreal &o) const {return x<=o.x;}
    bool operator > (const surreal &o) const {return x>o.x;}
    bool operator >= (const surreal &o) const {return x>=o.x;}

    surreal operator += (const surreal &o){
        x+=o.x;
        return *this;
    }
    surreal operator -= (const surreal &o){
        x-=o.x;
        return *this;
    }
    surreal operator + (const surreal &o){return surreal(*this)+=o;}
    surreal operator - (const surreal &o){return surreal(*this)-=o;}

    friend surreal avg(const surreal &a, const surreal &b){
        return surreal((a.x+b.x)/2);
    }

    friend ostream& operator << (ostream& os, const surreal &b){
        return os << b.x;
    }
};

template<typename T>
surreal<T> combine(surreal<T> l, surreal<T> r){
    assert(l<r);
    surreal<T> curl(-surreal<T>::inf),curr(surreal<T>::inf),mid(0);
    while(1){
        if(l<mid&&mid<r) return mid;
        if(r<=mid){
            curr=mid;
            if(curl<surreal<T>(-surreal<T>::inf/2)) mid=curr-1;
            else mid=avg(curl,curr);
        }
        else{
            curl=mid;
            if(curr>surreal<T>(surreal<T>::inf/2)) mid=curl+1;
            else mid=avg(curl,curr);
        }
    }
}
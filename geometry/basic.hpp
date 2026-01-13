#pragma once

template<typename T>
struct Point{
    inline static T eps;

    static void set_eps(T _eps){
        eps=_eps;
    }

    T x,y;
    Point():x(0),y(0){}

    template<typename T1, typename T2>
    Point(T1 _x, T2 _y):x(_x),y(_y){}

    template<typename T1, typename T2>
    Point(pair<T1,T2> p):x(p.first),y(p.second){}

    Point operator += (const Point &o){
        x+=o.x,y+=o.y;
        return *this;
    }
    Point operator -= (const Point &o){
        x-=o.x,y-=o.y;
        return *this;
    }

    Point operator + (const Point &o) const {return Point(x+o.x,y+o.y);}
    Point operator - (const Point &o) const {return Point(x-o.x,y-o.y);}
    Point operator - () const {return Point(-x,-y);}
    Point operator * (T k) {return Point(x*k,y*k);}
    Point operator / (T k) {return Point(x/k,y/k);}
    T operator * (Point o) {return x*o.x+y*o.y;}
    T operator ^ (Point o) {return x*o.y-y*o.x;}
    bool operator < (const Point &o) const {return x==o.x?y<o.y:x<o.x;}

    double angle(){return atan2(y,x);}

    friend istream& operator >> (istream& is, Point &b){
        T x,y;
        is >> x >> y;
        b=Point(x,y);
        return is;
    }
    friend ostream& operator << (ostream& os, const Point &b){
        return os << b.x << ' ' << b.y;
    }
};

template<typename T>
T abs2(Point<T> p){return p*p;}

template<typename T>
T abs(Point<T> p){
    return sqrt(abs2<T>(p));
}
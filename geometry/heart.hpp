#pragma once

// circumcenter
template<typename T>
Point<T> circenter(Point<T> p0, Point<T> p1, Point<T> p2){
    p1-=p0,p2-=p0;
    T x1=p1.x,y1=p1.y,x2=p2.x,y2=p2.y;   
    T m=T(2)*(x1*y2-y1*x2);
    Point<T> center;
    center.x=(x1*x1*y2-x2*x2*y1+y1*y2*(y1-y2))/m;
    center.y=(x1*x2*(x2-x1)-y1*y1*x2+x1*y2*y2)/m;
    return center+p0;
}
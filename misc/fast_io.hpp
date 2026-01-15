#pragma once

struct fast_io{
    static const int M=1<<15;
    char buf[M],*p=buf,*q=buf;
    char readchar(){
        if(p==q&&(q=(p=buf)+fread(buf,1,M,stdin))==buf) return EOF;
        else return *p++;
    }
    bool readint(int &x){
        int c=readchar(),neg=0;
        x=0;
        while((c<'0'||c>'9')&&c!='-'&&c!=EOF) c=readchar();
        if(c==EOF) return false;
        if(c=='-') neg=1,c=readchar();
        while(c>='0'&&c<='9') x=(x<<3)+(x<<1)+(c^'0'),c=readchar();
        if(neg) x=-x;
        return true;
    }
    bool readll(ll &x){
        int c=readchar(),neg=0;
        x=0;
        while((c<'0'||c>'9')&&c!='-'&&c!=EOF) c=readchar();
        if(c==EOF) return false;
        if(c=='-') neg=1,c=readchar();
        while(c>='0'&&c<='9') x=(x<<3)+(x<<1)+(c^'0'),c=readchar();
        if(neg) x=-x;
        return true;
    }
    bool readstr(string &a){
        a.clear();
        int c=readchar();
        while(c==' '||c=='\n') c=readchar();
        if(c==EOF) return false;
        while(c!=' '&&c!='\n'&&c!=EOF) a+=c,c=readchar();
        return true;
    }
    char buf2[M],*r=buf2;
    void writechar(char c){
        (*r++)=c;
        if(r==buf2+M) r=buf2,fwrite(buf2,1,M,stdout);
    }
    void writestr(string a, char end='\n'){
        for(char c: a) writechar(c);
        if(end) writechar(end);
    }
    void writeint(int x, char end='\n'){
        if(x<0) writechar('-'),x*=-1;
        if(x==0){
            writechar('0');
            return;
        }
        char tmp[14],*ptr=tmp;
        while(x){
            (*ptr++)='0'+(x%10);
            x/=10;
        }
        while(ptr!=tmp) writechar(*--ptr);
        if(end) writechar(end);
    }
    void writell(ll x, char end='\n'){
        if(x<0) writechar('-'),x*=-1;
        if(x==0){
            writechar('0');
            return;
        }
        char tmp[25],*ptr=tmp;
        while(x){
            (*ptr++)='0'+(x%10);
            x/=10;
        }
        while(ptr!=tmp) writechar(*--ptr);
        if(end) writechar(end);
    }
    void endl(){
        writechar('\n');
    }
    void end(){
        fwrite(buf2,1,r-buf2,stdout);
    }
}io;
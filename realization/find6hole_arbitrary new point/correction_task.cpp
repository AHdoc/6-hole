#include<iostream>
#include<cstdio>
#include<cstring>
#include<cmath>
#include<ctime>
#include<random>
#include<vector>
#include<set>
#include<queue>
#include<algorithm>
#include<assert.h>
using namespace std;

typedef long long LL;
LL absLL(LL x){return (x<0?-x:x);}
LL gcd(LL a,LL b){return (b==0?a:gcd(b,a%b));}

namespace ahdoc{
	#define x first
	#define y second
	#define c1 first.first
	#define c2 first.second
	#define c3 second
	#define MP3(a,b,c) make_pair(make_pair(a,b),c)
	
	const int MAXN=30;
	
	int n;
	pair<LL,LL> pt[MAXN];
	int Q_head[MAXN],Q_tail[MAXN];
	pair<pair<int,int>,int> Q[MAXN][MAXN];
	pair<pair<int,int>,int> C[MAXN];
	
	bool on_the_left(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return (b.x-a.x)*(c.y-a.y)-(c.x-a.x)*(b.y-a.y)>0;}
	
	bool proceed(int i,int j){
		//cerr<<">   i="<<i<<"   j="<<j<<"\n";
		if(on_the_left(pt[i],pt[j],make_pair(0LL,0LL))){
			while(Q_head[i]!=Q_tail[i] && on_the_left(pt[Q[i][Q_head[i]].c1],pt[i],pt[j])){
				if(proceed(Q[i][Q_head[i]].c1,j)) return true;
				auto tmp=Q[i][Q_head[i]];
				if(tmp.c1>C[i].c1) C[i].c1=tmp.c1;
				if(tmp.c2>C[i].c2) C[i].c2=tmp.c2;
				if(tmp.c3>C[i].c3) C[i].c3=tmp.c3;
				++Q_head[i];
			}
			if(C[i].c3!=-1 && on_the_left(pt[C[i].c3],pt[j],make_pair(0LL,0LL))) return true;
			Q[j][Q_tail[j]++]=MP3(i,C[i].c1,C[i].c2);
			//cerr<<"<   i="<<i<<"   j="<<j<<"\n";
			//for(int k=0;k<pt.size();k++){
			//	cerr<<"       Q["<<k<<"]={ ";
			//	for(auto x:Q[k]) cerr<<"("<<x.c1<<","<<x.c2<<","<<x.c3<<") ";
			//	cerr<<"}\n";
			//	cerr<<"       C["<<k<<"]=("<<C[k].c1<<","<<C[k].c2<<","<<C[k].c3<<")\n";
			//}
		}
		return false;
	}
	
	bool solve(){
		//cerr<<"n="<<n<<"\n";
		//for(int i=0;i<n;i++) cerr<<"   "<<i<<": ("<<pt[i].x<<","<<pt[i].y<<")\n";
		for(int i=0;i<n;i++){
			Q_head[i]=Q_tail[i]=0;
			C[i]=MP3(-1,-1,-1);
		}
		for(int i=0;i+1<n;i++)
			if(proceed(i,i+1))
				return true;
		return false; 
	}
	
	bool find6hole(vector<pair<LL,LL>> _pt,pair<LL,LL> p){
		n=_pt.size();
		for(int i=0;i<n;i++){
			pt[i].x=_pt[i].x-p.x;
			pt[i].y=_pt[i].y-p.y;
		}
		sort(pt,pt+n,[](pair<LL,LL> p1,pair<LL,LL> p2){return atan2(p1.y,p1.x)<atan2(p2.y,p2.x);});
		if(solve()) return true;
		for(int i=0;i<n;i++){
			pt[i].x=p.x-_pt[i].x;
			pt[i].y=p.y-_pt[i].y;
		}
		sort(pt,pt+n,[](pair<LL,LL> p1,pair<LL,LL> p2){return atan2(p1.y,p1.x)<atan2(p2.y,p2.x);});
		if(solve()) return true;
		return false;
	}
}

bool check(vector<pair<LL,LL>> pt,pair<LL,LL> p){
	set<pair<LL,LL>> ang;
	for(int j=0;j<pt.size();j++){
		LL x=pt[j].x-p.x,y=pt[j].y-p.y;
		if(x==0 && y==0) return false; // Two points coincide.
		if(y==0) return false; // Two points lie on a horizontal line.
		LL g=gcd(absLL(x),absLL(y)); x/=g; y/=g;
		if(x<0){x=-x; y=-y;}
		else if(x==0 && y<0) y=-y;
		ang.insert({x,y}); 
	}
	if(ang.size()!=pt.size()) return false;
	return !ahdoc::find6hole(pt,p);
}

LL crossproduct(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return (b.x-a.x)*(c.y-a.y)-(c.x-a.x)*(b.y-a.y);}
LL innerproduct(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return (b.x-a.x)*(c.x-a.x)+(b.y-a.y)*(c.y-a.y);}
LL crossproduct(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c,pair<LL,LL> d){return (b.x-a.x)*(d.y-c.y)-(d.x-c.x)*(b.y-a.y);}
bool on_the_left(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return crossproduct(a,b,c)>0;}
bool on_the_line(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return crossproduct(a,b,c)==0;}
void check_no6hole(vector<pair<LL,LL>> pt){ // Check whether or not the given set of points contains no 6-hole
	int n=pt.size();
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++){
			set<int> pts={i,j};
			if(pts.size()!=2) continue;
			if(pt[i].x==pt[j].x && pt[i].y==pt[j].y){
				cerr<<"Two points coincide.\n";
				exit(1);
			}
		}
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			for(int k=0;k<n;k++){
				set<int> pts={i,j,k};
				if(pts.size()!=3) continue;
				if(on_the_line(pt[i],pt[j],pt[k])){
					cerr<<"Three points which are colinear.\n";
					exit(1);
				}
			}
	int a[6];
	for(a[0]=0;a[0]<n;a[0]++) for(a[1]=0;a[1]<n;a[1]++) for(a[2]=0;a[2]<n;a[2]++) for(a[3]=0;a[3]<n;a[3]++) for(a[4]=0;a[4]<n;a[4]++) for(a[5]=0;a[5]<n;a[5]++){
		set<int> pts={a[0],a[1],a[2],a[3],a[4],a[5]};
		if(pts.size()!=6) continue;
		bool is6hole=true;
		for(int i=0;i<6;i++){
			int j=(i+1)%6;
			for(int k=(j+1)%6;k!=i;k=(k+1)%6)
				if(!on_the_left(pt[a[i]],pt[a[j]],pt[a[k]])){ is6hole=false; break; }
			if(!is6hole) break;
		}
		for(int k=0;k<n;k++){
			bool chk=true;
			for(int i=0;i<6;i++) if(a[i]==k){ chk=false; break; }
			if(!chk) continue;
			bool chkout=false;
			for(int i=0;i<6;i++){
				int j=(i+1)%6;
				if(!on_the_left(pt[a[i]],pt[a[j]],pt[k])){ chkout=true; break; }
			}
			if(!chkout){ is6hole=false; break; }
		}
		if(is6hole){
			cerr<<"n="<<n<<"\n";
			for(int k=0;k<n;k++) cerr<<"   Point #"<<k+1<<": ("<<pt[k].x<<","<<pt[k].y<<")\n";
			cerr<<"A 6-hole:\n";
			for(int i=0;i<6;i++) cerr<<"   "<<a[i]<<": ("<<pt[a[i]].x<<","<<pt[a[i]].y<<")\n";
			exit(1);
		}
	}
}

int main(){
	vector<pair<LL,LL>> pt;
	
	vector<pair<LL,LL>> table{{0,0},{397,302},{-649,-975},{-1072,470},{-540,-1845},{-130,1041},
	                          {524,-2017},{177,-1607},{-2086,2903},{3889,-1576},{-3289,-1786},
							  {2279,2948},{-254,3785},{1280,-2979},{-6681,-4707},{271,-346}};
	//sort(table.begin(),table.end());
	for(int i=0;i<16;i++){
		if(check(pt,table[i])) pt.push_back(table[i]);
		check_no6hole(pt);
	}
}


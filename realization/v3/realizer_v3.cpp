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

namespace ahdoc{
	typedef long long LL;
	
	#define x first
	#define y second
	
	const int MAXN=30;
	
	int n;
	pair<LL,LL> pt[MAXN];
	int Q_head[MAXN],Q_tail[MAXN];
	pair<int,int> Q[MAXN][MAXN];
	int C[MAXN];
	
	bool on_the_left(int a,int b,int c){return (pt[b].x-pt[a].x)*(pt[c].y-pt[a].y)-(pt[c].x-pt[a].x)*(pt[b].y-pt[a].y)>0;}
	
	bool proceed(int i,int j){
		while(Q_head[i]!=Q_tail[i] && on_the_left(Q[i][Q_head[i]].first,i,j)){
			if(proceed(Q[i][Q_head[i]].first,j)) return true;
			if(Q[i][Q_head[i]].second>C[i]) C[i]=Q[i][Q_head[i]].second;
			++Q_head[i];
		}
		if(C[i]>=3) return true;
		Q[j][Q_tail[j]++]=make_pair(i,C[i]+1);
		return false;
	}
	
	bool solve(){
		for(int i=0;i<n;i++){
			Q_head[i]=Q_tail[i]=0;
			C[i]=0;
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
		sort(pt,pt+n,[](pair<LL,LL> p1,pair<LL,LL> p2){return p1.x*p2.y-p1.y*p2.x>0;});
		if(solve()) return true;
		return false;
	}
}

namespace geo {
	typedef long double ld;
	
	const double eps=1e-3;
	const double maxr=1e9;
	
	default_random_engine gen;
	
	class P{
		public:
			double x,y;
			P(double x=0,double y=0):x(x),y(y){}
			friend ld Dot(P A,P B) {return A.x*B.x+A.y*B.y; }
			friend ld Length(P A) {return sqrt(Dot(A,A)); }
			friend P operator- (P A,P B) { return P(A.x-B.x,A.y-B.y); }
			friend P operator+ (P A,P B) { return P(A.x+B.x,A.y+B.y); }
			friend P operator* (P A,double p) { return P(A.x*p,A.y*p); }
	};
	
	typedef P V;
	double Cross(V A,V B) {return A.x*B.y - A.y*B.x;}
	
	struct Line{
		P p;
		V v;
		double ang;
		Line(){}
		Line(P p,V v):p(p),v(v) {ang=atan2(v.y,v.x); }
		bool operator<(const Line & L) const {
			return ang<L.ang;
		}
	};
	
	bool OnLeft(Line L,P p) {
		return Cross(L.v,p-L.p)>eps*Length(L.v);
	} 
	P GetIntersection(Line a,Line b) {
		V u=a.p-b.p;
		double t = Cross(b.v,u) / Cross(a.v, b.v);
		return a.p + a.v*t;
	}
	vector<P> HalfplaneIntersection(vector<Line> L) {
		sort(L.begin(),L.end());
		int n=L.size(),fi,la;
		vector<P> p; p.resize(n);
		vector<Line> q; q.resize(n);
		q[fi=la=0]=L[0];
		for(int i=1;i<=n-1;i++){
			while(fi<la && !OnLeft(L[i],p[la-1])) --la;
			while(fi<la && !OnLeft(L[i],p[fi])) ++fi;
			q[++la]=L[i];
			if(fabs(Cross(q[la].v, q[la-1].v))*maxr<eps*Length(q[la].v)*Length(q[la-1].v)){
				--la;
				if(OnLeft(q[la],L[i].p)) q[la]=L[i];
			}
			if (fi<la) p[la-1]=GetIntersection(q[la-1],q[la]);
		}
		while(fi<la && !OnLeft(q[fi],p[la-1])) --la;
		vector<P> res;
		if(la-fi<=1);
		else{
			p[la]=GetIntersection(q[la],q[fi]);
			for(int i=fi;i<=la;i++) res.push_back(p[i]);
		}
		return res;
	} 
	
	double PolygonArea(vector<P> p){
		double area=0.0;
		int n=p.size();
		for(int i=1;i+1<n;i++)
			area+=0.5*Cross(p[i]-p[0],p[i+1]-p[0]);
		return area;
	}
	double Realnum_Inside01(){
		uniform_real_distribution<> dist(0,1);
		return dist(gen);
	}
	int Weighted_Random(vector<double> weights){
		discrete_distribution<int> dist(weights.begin(),weights.end());
		return dist(gen);
	}
	P PointPicking_Triangle(vector<P> tri){
		double a=Realnum_Inside01(),b=Realnum_Inside01();
		if(a+b<1)
			return {tri[0].x+a*(tri[1].x-tri[0].x)+b*(tri[2].x-tri[0].x),tri[0].y+a*(tri[1].y-tri[0].y)+b*(tri[2].y-tri[0].y)};
		else
			return PointPicking_Triangle(tri);
	}
	P PointPicking_Polygon(vector<P> p){
		int n=p.size();
		vector<double> areas;
		for(int i=1;i+1<n;i++)
			areas.push_back(0.5*Cross(p[i]-p[0],p[i+1]-p[0]));
		int j=1+Weighted_Random(areas);
		return PointPicking_Triangle({p[0],p[j],p[j+1]}); 
	}
}

typedef long long LL;

#define x first
#define y second

LL absLL(LL x){return (x<0?-x:x);}
LL gcd(LL a,LL b){return (b==0?a:gcd(b,a%b));}

const int MAXN=30;

int n;
LL tot_query_check,tot_query_find6hole,achievement[MAXN+1];
LL radius[MAXN+1],lvl[MAXN+2];

LL crossproduct(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return (b.x-a.x)*(c.y-a.y)-(c.x-a.x)*(b.y-a.y);}
LL innerproduct(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return (b.x-a.x)*(c.x-a.x)+(b.y-a.y)*(c.y-a.y);}
LL crossproduct(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c,pair<LL,LL> d){return (b.x-a.x)*(d.y-c.y)-(d.x-c.x)*(b.y-a.y);}
bool on_the_left(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return crossproduct(a,b,c)>0;}
bool on_the_line(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return crossproduct(a,b,c)==0;}

long long tot_b=0; 
bool find6hole(vector<pair<LL,LL>> pt,pair<LL,LL> p){
	double ret=ahdoc::find6hole(pt,p);
	++tot_query_find6hole;
	tot_b+=ret;
	if(tot_query_find6hole%1000000LL==0){
		time_t current_t=time(NULL);
		time(&current_t);
		cerr<<"time:"<<ctime(&current_t)<<endl;
		cerr<<"tot_6holes:"<<tot_b<<endl;
		cerr<<"tot_query_find6hole = "<<tot_query_find6hole<<" & tot_query_check = "<<tot_query_check<<":";
		for(int j=1;j<=30;j++){
			if(achievement[j]==0) break;
			if(lvl[j-1]!=lvl[j]) cerr<<"\n   lvl="<<lvl[j]<<"   ";
			cerr<<j<<":"<<achievement[j]<<" ";
		}
		cerr<<"\n";
	}
	return ret;
}

vector<geo::P> get_range(int i,int ii,vector<pair<LL,LL>> pt){
	geo::P A((double)radius[i],(double)radius[i]),B((double)-radius[i],(double)radius[i]),C((double)-radius[i],(double)-radius[i]),D((double)radius[i],(double)-radius[i]);
	vector<geo::Line> L={geo::Line(A,B-A),geo::Line(B,C-B),geo::Line(C,D-C),geo::Line(D,A-D)};
	if(lvl[i]==lvl[i-1]+1){
		for(int j=1;j<i;j++)
			L.push_back(geo::Line(geo::P((double)pt[j-1].x+1,0.0),geo::P(0.0,-1.0)));
	}else if(lvl[i]+1==lvl[i+1]){
		for(int j=1;j<i;j++){
			if(j!=i-1)
				L.push_back(geo::Line(geo::P(pt[j-1].x,pt[j-1].y),geo::P(pt[i-2].x,pt[i-2].y)-geo::P(pt[j-1].x,pt[j-1].y)));
			if(j!=ii)
				L.push_back(geo::Line(geo::P(pt[ii-1].x,pt[ii-1].y),geo::P(pt[j-1].x,pt[j-1].y)-geo::P(pt[ii-1].x,pt[ii-1].y)));
		}
	}else{
		for(int j=1;j<i;j++)
			if(j!=i-1)
				L.push_back(geo::Line(geo::P(pt[j-1].x,pt[j-1].y),geo::P(pt[i-2].x,pt[i-2].y)-geo::P(pt[j-1].x,pt[j-1].y)));
	}
	return geo::HalfplaneIntersection(L);
}

bool check(int i,int ii,vector<pair<LL,LL>> pt,pair<LL,LL> p){
	++tot_query_check;
	set<pair<LL,LL>> ang;
	for(int j=1;j<i;j++){
		LL x=pt[j-1].x-p.x,y=pt[j-1].y-p.y;
		if(x==0 || y==0) return false;
		LL g=gcd(absLL(x),absLL(y)); x/=g; y/=g;
		if(x<0){x=-x; y=-y;}
		else if(x==0) y=-y;
		ang.insert({x,y}); 
	}
	if(ang.size()!=i-1) return false;
	if(lvl[i]==lvl[i-1]+1){
		for(int j=1;j<i;j++)
			if(p.x<=pt[j-1].x) return false;
	}else if(lvl[i]+1==lvl[i+1]){
		for(int j=1;j<i;j++)
			if((j!=i-1 && !on_the_left(p,pt[j-1],pt[i-2])) || (j!=ii && !on_the_left(p,pt[ii-1],pt[j-1]))) return false;
	}else{
		for(int j=1;j<i;j++)
			if(j!=i-1 && !on_the_left(p,pt[j-1],pt[i-2])) return false;
	}
	return !find6hole(pt,p);
}

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

pair<LL,LL> random_walk(int i,pair<LL,LL> p){
	vector<pair<LL,LL>> dir={{0,1},{0,-1},{1,0},{-1,0}};
	int o=rand()%4;
	p.x+=dir[o].x; p.y+=dir[o].y;
	return p;
}

int depthdiff_to_confidence(int i,int x){
	if(x>=1) return x;
	else return 1;
}

LL dfs(int i,int ii,vector<pair<LL,LL>> pt){ // points numbered from ii to i are of the same level
	++achievement[i-1];
	if(i==n+1){
		for(int i=1;i<=n;i++) cerr<<"Point #"<<i<<": ("<<pt[i-1].x<<","<<pt[i-1].y<<")\n";
		check_no6hole(pt);
		exit(0);
	}
	
	LL max_depth=i;
	const LL initlvl=2;
	const LL base=2;
	LL amo=1;
	switch((lvl[i]-(lvl[i-1]!=lvl[i]))*(lvl[i-1]!=lvl[i] || lvl[i]!=lvl[i+1])){
		case initlvl: amo=base; break;
		case initlvl+1: amo=base*base; break;
		case initlvl+2: amo=base*base*base; break;
		case initlvl+3: amo=base*base*base*base; break;
		case initlvl+4: amo=base*base*base*base*base; break;
		case initlvl+5: amo=base*base*base*base*base*base; break;
		case initlvl+6: amo=base*base*base*base*base*base; break;
		case initlvl+7: amo=base*base*base*base*base*base; break;
		case initlvl+8: amo=base*base*base*base*base*base; break;
		case initlvl+9: amo=base*base*base*base*base*base; break;
	}
	vector<geo::P> range=get_range(i,ii,pt); 
	if(geo::PolygonArea(range)<geo::eps) amo=0;
	for(LL t=1;t<=amo;t++){
		geo::P _p=geo::PointPicking_Polygon(range); 
		pair<LL,LL> p=make_pair((LL)_p.x,(LL)_p.y);
		LL max_depthdiff=0;
		for(LL t_conf=1;t_conf<=depthdiff_to_confidence(i,max_depthdiff);t_conf++){
			if(t_conf==1){
				if(!check(i,ii,pt,p)) continue;
			}else{
				pair<LL,LL> p2=random_walk(i,p);
				if(!check(i,ii,pt,p2)) continue;
				else p=p2;
			}
			vector<pair<LL,LL>> pt2=pt;
			pt2.push_back(p);
			//check_no6hole(pt2); // bruteforce
			max_depthdiff=max(max_depthdiff,dfs(i+1,lvl[i]==lvl[i+1]?ii:i+1,pt2)-i);
			max_depth=max(max_depth,i+max_depthdiff);
		}
	}
	return max_depth;
}

void Realizer(string pat){
	n=0;
	LL radius0=1LL;
	for(int i=pat.size()-1,j=1,k;i>=0;i--,j=k){
		k=j+pat[i]-'0';
		for(int l=j;l<k;l++){
			lvl[l]=pat.size()-1-i;
			radius[l]=radius0;
		}
		radius0*=10LL;
		n=k-1;
	}
	lvl[0]=0; lvl[n+1]=pat.size();
	cerr<<"pat="<<pat<<"   n="<<n<<"\n"; 
	for(int i=1;i<=n;i++) cerr<<"i="<<i<<":   radius="<<radius[i]<<"   lvl="<<lvl[i]<<"\n";
	cerr<<"lvl[0]="<<lvl[0]<<"   lvl[n+1]="<<lvl[n+1]<<"\n";
	
	tot_query_find6hole=0; for(int i=1;i<=30;i++) achievement[i]=0;
	for(;;) dfs(2,lvl[1]==lvl[2]?1:2,{{0LL,0LL}}); // must fix the kernel point
}

int main(){
	//Realizer("558620");
	//Realizer("346650");
	
	//Realizer("333330"); //done with base=2, 1s.
	//Realizer("3333330"); //done with base=2, 8s. 
	Realizer("4444440");
	
	//Realizer("8730");
	//Realizer("88510");
	//Realizer("3477710");
	
	//check({{0,0},{59,-35},{-99,81},{-77,6},{16,-87},{96,-82}});
	//cout<<ahdoc::find6hole({{0,0},{59,-35},{-99,81},{-77,6},{16,-87},{96,-82}},{92,-73})<<"\n";
	//check({{0,0},{59,-35},{-99,81},{-77,6},{16,-87},{96,-82},{92,-73}});
}


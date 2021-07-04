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
	#define c1 first.first
	#define c2 first.second
	#define c3 second
	#define MP3(a,b,c) make_pair(make_pair(a,b),c)
	
	bool on_the_left(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return (b.x-a.x)*(c.y-a.y)-(c.x-a.x)*(b.y-a.y)>0;}
	
	vector<pair<LL,LL>> pt;
	vector<queue<pair<pair<int,int>,int>>> Q;
	vector<pair<pair<int,int>,int>> C;
	
	bool proceed(int i,int j){
		//cerr<<">   i="<<i<<"   j="<<j<<"\n";
		while(!Q[i].empty() && on_the_left(pt[Q[i].front().c1],pt[i],pt[j])){
			if(proceed(Q[i].front().c1,j)) return true;
			auto tmp=Q[i].front();
			if(tmp.c1>C[i].c1) C[i].c1=tmp.c1;
			if(tmp.c2>C[i].c2) C[i].c2=tmp.c2;
			if(tmp.c3>C[i].c3) C[i].c3=tmp.c3;
			Q[i].pop();
		}
		if(C[i].c3!=-1 && on_the_left(pt[C[i].c3],pt[j],make_pair(0LL,0LL))) return true;
		Q[j].push(MP3(i,C[i].c1,C[i].c2));
		//cerr<<"<   i="<<i<<"   j="<<j<<"\n";
		//for(int k=0;k<pt.size();k++){
		//	cerr<<"       Q["<<k<<"]={ ";
		//	for(auto x:Q[k]) cerr<<"("<<x.c1<<","<<x.c2<<","<<x.c3<<") ";
		//	cerr<<"}\n";
		//	cerr<<"       C["<<k<<"]=("<<C[k].c1<<","<<C[k].c2<<","<<C[k].c3<<")\n";
		//}
		return false;
	}
	
	bool solve(vector<pair<LL,LL>> _pt){
		pt=_pt;
		int n=pt.size();
		//cerr<<"n="<<n<<"\n";
		//for(int i=0;i<n;i++) cerr<<"   "<<i<<": ("<<pt[i].x<<","<<pt[i].y<<")\n";
		Q.clear(); C.clear(); Q.resize(n); C.resize(n);
		for(int i=0;i<n;i++)
			C[i].c1=C[i].c2=C[i].c3=-1;
		for(int i=0;i+1<n;i++)
			if(proceed(i,i+1))
				return true;
		return false; 
	}
	
	bool find6hole(vector<pair<LL,LL>> pt,pair<LL,LL> p){
		int n=pt.size();if (n<5) return 0;
		for(int i=0;i<n;i++){
			pt[i].x-=p.x;
			pt[i].y-=p.y;
		}
		sort(pt.begin(),pt.end(),[](pair<LL,LL> p1,pair<LL,LL> p2){return atan2(p1.y,p1.x)<atan2(p2.y,p2.x);});
		if(solve(pt)) return true;
		for(int i=0;i<n;i++){
			pt[i].x*=-1;
			pt[i].y*=-1;
		}
		sort(pt.begin(),pt.end(),[](pair<LL,LL> p1,pair<LL,LL> p2){return atan2(p1.y,p1.x)<atan2(p2.y,p2.x);});
		if(solve(pt)) return true;
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
	double Realnum_Inside01(){return (double)(rand()+1)/32769;}
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
LL random(LL a,LL b){
	assert(a<=b);
	LL x=rand()*32768;
	x=(x+rand())*32768;
	x=(x+rand())*32768;
	x=(x+rand())%(b-a+1);
	x+=a;
	return x;
}

const int MAXN=30;

int n;
LL tot_query_check,tot_query_find6hole,achievement[MAXN+1];
LL radius[MAXN+1],lvl[MAXN+2];

LL crossproduct(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return (b.x-a.x)*(c.y-a.y)-(c.x-a.x)*(b.y-a.y);}
LL innerproduct(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return (b.x-a.x)*(c.x-a.x)+(b.y-a.y)*(c.y-a.y);}
LL crossproduct(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c,pair<LL,LL> d){return (b.x-a.x)*(d.y-c.y)-(d.x-c.x)*(b.y-a.y);}
bool on_the_left(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return crossproduct(a,b,c)>0;}
bool on_the_line(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return crossproduct(a,b,c)==0;}
namespace geo_ll{
	#define For(i,n) for(int i=1;i<=n;i++)
	#define Fork(i,k,n) for(int i=k;i<=n;i++)
	#define Rep(i,n) for(int i=0;i<n;i++)
	#define ForD(i,n) for(int i=n;i;i--)
	#define ForkD(i,k,n) for(int i=n;i>=k;i--)
	#define RepD(i,n) for(int i=n;i>=0;i--)
	#define Forp(x) for(int p=Pre[x];p;p=Next[p])
	#define Forpiter(x) for(int &p=iter[x];p;p=Next[p])  
	#define Lson (o<<1)
	#define Rson ((o<<1)+1)
	#define MEM(a) memset(a,0,sizeof(a));
	#define MEMI(a) memset(a,127,sizeof(a));
	#define MEMi(a) memset(a,128,sizeof(a));
	#define INF (2139062143)
	#define F (100000007)
	#define pb push_back
	#define mp make_pair 
	#define fi first
	#define se second
	#define vi vector<int> 
	#define pi pair<int,int>
	#define SI(a) ((a).size())
	#define ALL(x) (x).begin(),(x).end()
	typedef long long ll;
	ll mul(ll a,ll b){return (a*b)%F;}
	ll add(ll a,ll b){return (a+b)%F;}
	ll sub(ll a,ll b){return (a-b+llabs(a-b)/F*F+F)%F;}
	void upd(ll &a,ll b){a=(a%F+b%F)%F;}
	int read()
	{
		int x=0,f=1; char ch=getchar();
		while(!isdigit(ch)) {if (ch=='-') f=-1; ch=getchar();}
		while(isdigit(ch)) { x=x*10+ch-'0'; ch=getchar();}
		return x*f;
	} 
	ll sqr(ll a){return a*a;}
	ll dcmp(ll x){
		if(x>0) return 1;if(x<0) return -1;return 0;
	}
	class P{
		public:
			ll x,y;
			P(ll x=0,ll y=0):x(x),y(y){}
			
			friend ll dis2(P A,P B){return sqr(A.x-B.x)+sqr(A.y-B.y);	}
			friend ll dis2(P A){return sqr(A.x)+sqr(A.y);	}
			friend ll Dot(P A,P B) {return A.x*B.x+A.y*B.y; }
			friend P operator- (P A,P B) { return P(A.x-B.x,A.y-B.y); }
			P(P A,P B):x(B.x-A.x),y(B.y-A.y){}
			friend P operator+ (P A,P B) { return P(A.x+B.x,A.y+B.y); }
			friend P operator* (P A,double p) { return P(A.x*p,A.y*p); }
			friend P operator/ (P A,double p) { return P(A.x/p,A.y/p); }
			friend bool operator< (const P& a,const P& b) {return dcmp(a.x-b.x)<0 ||(dcmp(a.x-b.x)==0&& dcmp(a.y-b.y)<0 );}
		}; 
	P read_point() {
		P a;
		scanf("%lld%lld",&a.x,&a.y);
		return a;	
	} 
	bool operator==(const P& a,const P& b) {
		return dcmp(a.x-b.x)==0 && dcmp(a.y-b.y) == 0;
	} 
	typedef P V;
	
	ll Cross(V A,V B) {return A.x*B.y - A.y*B.x;}
	ll Area2(P A,P B,P C) {return Cross(B-A,C-A);}
	
	bool OnLeft(P A,P B,P C) {
		return Cross(B-A,C-A)>0;
	} 
	int Quadrant(P a)
	{
	    if(a.x>0&&a.y>=0) return 1;
	    if(a.x<=0&&a.y>0) return 2;
	    if(a.x<0&&a.y<=0) return 3;
	    return 4;
	}
	P _p;
	int cmp(P A,P B) //1:a>b 0:a<=b
	{
		if(Quadrant(A-_p)!=Quadrant(B-_p))
			return Quadrant(A-_p)<Quadrant(B-_p);
		ll tmp=Cross(V(_p,A),V(_p,B));
		if (tmp>0) return 1;
		else if (tmp==0) return (-(dis2(_p,A)-dis2(_p,B))>0)?1:0;
		else return 0;
	}
	void PolarSort(vector<P> &v,P _p2=P(0,0)) {
		_p=_p2;
		sort(ALL(v),cmp);
	}

	int cmp_O(P A,P B) //1:a>b 0:a<=b
	{
		
		if(Quadrant(A)!=Quadrant(B))
			return Quadrant(A)<Quadrant(B);
		ll tmp=Cross(A,B);
		if (tmp>0) return 1;
		else if (tmp==0) return (-(dis2(A)-dis2(B))>0)?1:0;
		else return 0;
	}
	
	struct Line{
		P p;
		V v;
		double ang;
		Line(){}
		Line(P p,V v):p(p),v(v) {ang=atan2(v.y,v.x); }
		bool operator<(const Line & L) const {
			return ang<L.ang;
		}
		P point(double a) {
			return p+v*a;
		}
	};
	bool OnLeft(Line L,P p) {
		return Cross(L.v,p-L.p)>0;
	} 
	bool OnRight(Line L,P p) {
		return Cross(L.v,p-L.p)<0;
	} 
	class Find6hole{
	public:
		int n;
		#define MAXN (31)

		pair<pair<int,int>, int>  q[MAXN][MAXN];
		pair<pair<int,int>, int>  C[MAXN];
		int q_h[MAXN],q_t[MAXN];
		P _p;
		P vp[MAXN];
		
		bool proceed(int i,int j) {
//			cerr<<i<<' '<<j<<endl;
			while(q_h[i]<=q_t[i] && OnLeft(Line(vp[q[i][q_h[i]].fi.fi],vp[i]-vp[q[i][q_h[i]].fi.fi]),vp[j])){
			 //if k can see i && i can see j && turn)left, then k can see j
				if (proceed(q[i][q_h[i]].fi.fi,j)) return 1; // add k-j and p-k-j 
					
					auto now=q[i][q_h[i]];
					C[i].fi.fi=max(C[i].fi.fi,now.fi.fi);
					C[i].fi.se=max(C[i].fi.se,now.fi.se);
					C[i].se=max(C[i].se,now.se);
					q_h[i]++;
			}
			if(C[i].se!=-1 && OnLeft(Line(vp[C[i].se],vp[j]-vp[C[i].se]),_p ) ) {
				return 1;
			} 
			q[j][++q_t[j]]=(mp(mp(i,C[i].fi.fi),C[i].fi.se));
			return 0;
		}
//		
//		bool proceed(int i,int j){
//			//cerr<<">   i="<<i<<"   j="<<j<<"\n";
//			while(!q[i].empty() && OnLeft(vp[q[i].front().fi.fi],vp[i],vp[j])){
//				if(proceed(q[i].front().fi.fi,j)) return true;
//				auto tmp=q[i].front();
//				if(tmp.c1>C[i].c1) C[i].c1=tmp.c1;
//				if(tmp.c2>C[i].c2) C[i].c2=tmp.c2;
//				if(tmp.c3>C[i].c3) C[i].c3=tmp.c3;
//				q[i].pop();
//			}
//			if(C[i].c3!=-1 && OnLeft(vp[C[i].c3],vp[j],P(0LL,0LL))) return true;
//			q[j].push(mp( mp(i,C[i].c1),C[i].c2));
//			//cerr<<"<   i="<<i<<"   j="<<j<<"\n";
//			//for(int k=0;k<pt.size();k++){
//			//	cerr<<"       Q["<<k<<"]={ ";
//			//	for(auto x:Q[k]) cerr<<"("<<x.c1<<","<<x.c2<<","<<x.c3<<") ";
//			//	cerr<<"}\n";
//			//	cerr<<"       C["<<k<<"]=("<<C[k].c1<<","<<C[k].c2<<","<<C[k].c3<<")\n";
//			//}
//			return false;
//		}
		bool find6hole_R(){
			//Ci_1,Ci_2,Ci_3 Ci_j:Then chain starts in C_i, and end in i, length j
//			cout<<n<<endl;
			Rep(i,n) C[i]=mp(mp(-1,-1),-1);
			Rep(i,n) q_h[i]=0,q_t[i]=-1;
			Rep(i,n-1) {
				if(proceed(i,i+1))return 1;
			}
			return 0;			
		}
//		bool find6hole(vector<P> _vp,P __p){
//			vp=_vp;
//			_p=__p;
//			n=vp.size();
//			if (n<5) return 0; 
//			for(int i=0;i<n;i++){
//				vp[i]=vp[i]-_p;
//			}
//			_p=P(0,0);
//			PolarSort(vp,_p);
//			if(find6hole_R()) return 1;
//			for(int i=0;i<n;i++) vp[i].x*=-1,vp[i].y*=-1;
//			PolarSort(vp,_p);
//			if(find6hole_R()) return 1;
//			return 0;
//		}
		bool find6hole(vector<P> _vp,P __p){
			_p=__p;
			n=_vp.size();
			if (n<5) return 0;
			for(int i=0;i<n;i++){
				vp[i]=_vp[i]-_p;
			}
			_p=P(0,0);
			sort(vp,vp+n,cmp_O);
//			PolarSort(vp,_p);
			if(find6hole_R()) return 1;
//			int ans=-1,l=1,r=n-1;
//			while(l<=r) {
//				int m=l+r>>1;
//				if(Quadrant(vp[m])>=3) ans=m,l=m+1;else r=m-1;
//			}
//			if(ans!=-1){
//				rotate(vp.begin(),vp.begin()+ans,vp.end());
//				return find6hole_R();
//			}
			return 0;
		}
	}S;
}
long long tot_b=0; 
bool find6hole(vector<pair<LL,LL>> pt,pair<LL,LL> p){
	vector<geo_ll::P> vp;
	for(auto p:pt) {
		vp.pb(geo_ll::P(p.x,p.y));
	}
	geo_ll::P _p=geo_ll::P(p.x,p.y);
	
	bool b=geo_ll::S.find6hole(vp,_p);
	
//	bool b=ahdoc::find6hole(pt,p);
	++tot_query_find6hole;
	tot_b+=b;
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
	return b;
	
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
	const LL base=10;
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
	for(LL t=1;t<=2*amo;t++){
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
//			check_no6hole(pt2);
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
	Realizer("346650");
	
	//Realizer("333330"); //done with base=2, <5s.
	//Realizer("3333330"); //done with base=2, <30s. 
	
	//Realizer("8730");
	//Realizer("88510");
	//Realizer("3477710");
	
	//check({{0,0},{59,-35},{-99,81},{-77,6},{16,-87},{96,-82}});
	//cout<<ahdoc::find6hole({{0,0},{59,-35},{-99,81},{-77,6},{16,-87},{96,-82}},{92,-73})<<"\n";
	//check({{0,0},{59,-35},{-99,81},{-77,6},{16,-87},{96,-82},{92,-73}});
}


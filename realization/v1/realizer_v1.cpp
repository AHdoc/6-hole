#include<iostream>
#include<cstdio>
#include<cstring>
#include<cmath>
#include<ctime>
#include<random>
#include<vector>
#include<set>
#include<algorithm>

using namespace std;

namespace geo {
	typedef long double ld;
	
	const double eps=1e-6;
	const double maxr=1e6;
	
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
}

/**************************************************************/

struct Tpoint{
	double x,y;
	Tpoint(){}
	Tpoint(double _x,double _y){x=_x; y=_y;}
	Tpoint operator -(const Tpoint &b)const{return Tpoint(x-b.x,y-b.y);}
	double operator ^(const Tpoint &b)const{return x*b.y-y*b.x;}
	double operator *(const Tpoint &b)const{return x*b.x+y*b.y;}
};
double norm(Tpoint p){return sqrt(p.x*p.x+p.y*p.y);}

struct Tline{
	Tpoint s,e;
	double k;
	Tline(){}
	Tline(Tpoint _s,Tpoint _e){
		s=_s; e=_e;
		k=atan2(e.y-s.y,e.x-s.x);
	}
	Tpoint operator &(const Tline &b)const{
		Tpoint res=s;
		double t=((s-b.s)^(b.s-b.e))/((s-e)^(b.s-b.e));
		res.x+=(e.x-s.x)*t;
		res.y+=(e.y-s.y)*t;
		return res;
	}
};

double eps=1e-6;

default_random_engine gen;

vector<Tpoint> HPI(vector<Tline> line){
	int n=line.size();
	vector<geo::Line> l;
	for(int i=0;i<n;i++)
		l.push_back(geo::Line(geo::P(line[i].s.x,line[i].s.y),geo::P(line[i].e.x-line[i].s.x,line[i].e.y-line[i].s.y)));
	vector<geo::P> p=geo::HalfplaneIntersection(l);
	
	vector<Tpoint> res;
	if(!p.empty()){
		res.push_back(Tpoint(p[0].x,p[0].y));
		for(int i=1;i<p.size();i++)
			if(norm(res[res.size()-1]-Tpoint(p[i].x,p[i].y))>eps) res.push_back(Tpoint(p[i].x,p[i].y));
		if(norm(res[0]-res[res.size()-1])<eps) res.pop_back();
		if(res.size()<3) res.clear();
	}
	
	double area=0.0,area_triangle;
	for(int i=2;i<res.size();i++){ // Triangle (0,i-1,i)
		area_triangle=((res[i-1]-res[0])^(res[i]-res[0]));
		if(area_triangle<0){
			cerr<<"A triangle of negative area "<<fixed<<area_triangle<<" appears in the polygon of an HPI problem.\n";
			cerr<<n<<"\nLines:\n";
			for(int j=0;j<n;j++) cerr<<"   "<<fixed<<line[j].s.x<<" "<<line[j].s.y<<" "<<line[j].e.x<<" "<<line[j].e.y<<"\n";
			cerr<<"Points:\n";
			for(int j=0;j<res.size();j++) cerr<<"   "<<fixed<<res[j].x<<" "<<res[j].y<<"\n";
			int j;
			j=0;	cerr<<fixed<<res[j].x<<" "<<res[j].y<<"\n";
			j=i-1;	cerr<<fixed<<res[j].x<<" "<<res[j].y<<"\n";
			j=i;	cerr<<fixed<<res[j].x<<" "<<res[j].y<<"\n";
			res.clear();
			exit(1);
		}
		area+=area_triangle;
	}
	//if(!res.empty() && area<eps) cerr<<"A super small polygon appears with the area="<<area<<".\n"; 
	return res;
}

double Realnum_Inside01(){return (double)(rand()+1)/32769;}

int Weighted_Random(vector<double> weights){
	discrete_distribution<int> dist(weights.begin(),weights.end());
	return dist(gen);
}

Tpoint PointPicking_Triangle(vector<Tpoint> tri){
	double a=Realnum_Inside01(),b=Realnum_Inside01();
	if(a+b<1)
		return {tri[0].x+a*(tri[1].x-tri[0].x)+b*(tri[2].x-tri[0].x),tri[0].y+a*(tri[1].y-tri[0].y)+b*(tri[2].y-tri[0].y)};
	else
		return PointPicking_Triangle(tri);
}

Tpoint PointPicking_Polygon(vector<Tpoint> pol){
	vector<double> areas;
	for(int i=2;i<pol.size();i++) areas.push_back((pol[i-1]-pol[0])^(pol[i]-pol[0])); // Triangle (0,i-1,i)
	int j=2+Weighted_Random(areas);
	return PointPicking_Triangle({pol[0],pol[j-1],pol[j]}); 
}

const int MAXN=30;
const double MAXR=1e6;

int n,tot_achievement,achievement[MAXN+1];
double radius[MAXN+1];
int id[MAXN+1],idx[MAXN+1][MAXN+1][MAXN+1];
bool var[MAXN*MAXN*MAXN+1];

bool query(int _i,int _j,int _k){
	int i=id[_i],j=id[_j],k=id[_k];
	if(idx[i][j][k]>0) return var[idx[i][j][k]];
	else return (!var[-idx[i][j][k]]);
}

void dfs(int i,vector<Tpoint> pt){
	++tot_achievement;
	++achievement[i-1];
	if(tot_achievement%1000000==0){
		cerr<<"achieve = "<<tot_achievement<<"     ";
		for(int j=10;j<=30;j++){
			if(achievement[j]==0) break;
			cerr<<j<<":"<<achievement[j]<<" ";
		}
		cerr<<"\n";
	}
	
	if(i==n+1){
		for(int i=1;i<=n;i++) cerr<<fixed<<"Point #"<<i<<": Tpoint("<<pt[i-1].x<<","<<pt[i-1].y<<")\n";
		// Check
		for(int i=1;i<=n;i++)
			for(int j=1;j<=n;j++)
				for(int k=1;k<=n;k++){
					set<int> ijk={i,j,k};
					if(ijk.size()!=3) continue;
					double crossprod=(pt[j-1]-pt[i-1])^(pt[k-1]-pt[i-1]);
					if(crossprod>0 && query(i,j,k));
					else if(crossprod<0 && (!query(i,j,k)));
					else{
						cerr<<fixed<<"point i = ("<<pt[i-1].x<<","<<pt[i-1].y<<")\n";
						cerr<<fixed<<"point j = ("<<pt[j-1].x<<","<<pt[j-1].y<<")\n";
						cerr<<fixed<<"point k = ("<<pt[k-1].x<<","<<pt[k-1].y<<")\n";
						cerr<<fixed<<"crossprod = "<<crossprod<<"\n";
						cerr<<"query(i,j,k) = "<<query(i,j,k)<<"\n";
						exit(1);
					}
				}
		exit(0);
	}
	
	Tpoint A(-radius[i],-radius[i]),B(radius[i],-radius[i]),C(radius[i],radius[i]),D(-radius[i],radius[i]);
	vector<Tline> line={Tline(A,B),Tline(B,C),Tline(C,D),Tline(D,A)};
	for(int j=1;j<i;j++)
		for(int k=j+1;k<i;k++)
			if(query(j,k,i)) line.push_back(Tline(pt[j-1],pt[k-1]));
			else line.push_back(Tline(pt[k-1],pt[j-1]));
	vector<Tpoint> res=HPI(line);
	if(res.size()<3) return;
	for(int t=1;t<=2;t++){
		vector<Tpoint> pt2=pt;
		pt2.push_back(PointPicking_Polygon(res));
		dfs(i+1,pt2);
	}
}

void Realizer(string pat){
	freopen((pat+".txt").c_str(),"r",stdin);
	cin>>n;
	for(int i=1;i<=n;i++) radius[i]=MAXR;
	//double radius0=MAXR;
	//for(int i=0,j=n,k;i<pat.size();i++,j=k){
	//	k=j-pat[i]-'0';
	//	for(int l=j;l>k;l--) radius[l]=radius0;
	//	radius0=sqrt(radius0);
	//}
	for(int i=1;i<=n;i++) id[i]=n+1-i;
	
	int tot=0;
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				int x=++tot;
				idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=x;
				idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-x;
			}
	int x;
	while(cin>>x && x!=0){
		if(abs(x)<=tot){
			if(x>0) var[x]=true;
			else var[abs(x)]=false;
		}
	}
	fclose(stdin);
	
	tot_achievement=0; for(int i=1;i<=30;i++) achievement[i]=0;
	Tpoint P1(0,0); 
	for(int i=1;;i++) dfs(2,{P1});
}

int main(){
	cerr.precision(17); cout.precision(17);
	//Realizer("8730");
	Realizer("3477710");
}


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <random>
#include <vector>
#include <set>
#include <string>
#include <memory>
#include <queue>
#include <algorithm>
#include <assert.h>
#include <omp.h>
using namespace std;

typedef long long LL;
#define x first
#define y second
typedef pair<LL, LL> Point;
const int MAXN = 30;

#include <chrono>
using namespace std::chrono;
LL gettime() {
	return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

LL randint(LL l, LL r) {
	thread_local mt19937_64 generator(omp_get_thread_num());
	// thread_local mt19937_64 generator(123);
	return uniform_int_distribution<LL>(l, r)(generator);
}

namespace ahdoc{
	typedef LL LL;
	#define x first
	#define y second
	#define c1 first.first
	#define c2 first.second
	#define c3 second
	#define MP3(a,b,c) make_pair(make_pair(a,b),c)
	
	const int MAXN=30;
	
	thread_local int n;
	thread_local Point pt[MAXN];
	thread_local int Q_head[MAXN],Q_tail[MAXN];
	thread_local pair<pair<int,int>,int> Q[MAXN][MAXN];
	thread_local pair<pair<int,int>,int> C[MAXN];

	bool atan2_less(const Point &p1, const Point &p2) {
		// should be the same as atan2(p1.y,p1.x)<atan2(p2.y,p2.x);
		// not working when any of {p1.x, p1.y, p2.x, p2.y} is 0 !!!
		//    2 |  1
		//s: ---o-->
		//   -2 | -1
		int s1 = (((p1.y >= 0) << 1) - 1) << (p1.x < 0);
		int s2 = (((p2.y >= 0) << 1) - 1) << (p2.x < 0);
		if (s1 != s2)
			return s1 < s2;
		return (p1.x * p2.y - p2.x * p1.y) > 0;
	}
	
	bool on_the_left(const Point &a,const Point &b,const Point &c){return (b.x-a.x)*(c.y-a.y)-(c.x-a.x)*(b.y-a.y)>0;}
	bool on_the_line(const Point &a,const Point &b,const Point &c){return (b.x-a.x)*(c.y-a.y)-(c.x-a.x)*(b.y-a.y)==0;}
	
	bool proceed(int i,int j){
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
		}
		return false;
	}
	
	bool solve(){
		for(int i=0;i<n;i++){
			Q_head[i]=Q_tail[i]=0;
			C[i]=MP3(-1,-1,-1);
		}
		for(int i=0;i+1<n;i++)
			if(proceed(i,i+1))
				return true;
		return false; 
	}
	
	bool check(const vector<Point> &pt_, const Point &p){
		n = pt_.size();
		for (int i = 0; i < n; i++) {
			pt[i].x = pt_[i].x - p.x;
			pt[i].y = pt_[i].y - p.y;
			if (pt[i].x == 0 && pt[i].y == 0)
				return false;
		}
		sort(pt, pt + n, atan2_less);
		for (int i = 1; i < n; i++)
			if (pt[i - 1].x * pt[i].y == pt[i].x * pt[i - 1].y)
				return false;
		if (solve()) return false;
		for (int i = 0; i < n; i++) {
			pt[i].x = -pt[i].x;
			pt[i].y = -pt[i].y;
		}
		sort(pt, pt + n, atan2_less);
		if (solve()) return false;
		return true;
	}

	bool check_no6hole(vector<Point> pt){ 
		// Check whether or not the given set of points contains no 6-hole
		sort(pt.begin(), pt.end());
		vector<Point> pt2 = {pt[0]};
		for(int i=1;i<pt.size();i++)
			if(check(pt2,pt[i]))
				pt2.push_back(pt[i]);
			else
				return false;
		return true;
	}

	bool check_no6hole_bruteforce(vector<Point> pt, bool veb=false){ // Check whether or not the given set of points contains no 6-hole
		int n=pt.size();
		for(int i=0;i<n;i++)
			for(int j=0;j<n;j++){
				set<int> pts={i,j};
				if(pts.size()!=2) continue;
				if(pt[i].x==pt[j].x && pt[i].y==pt[j].y){
					if (veb)std::cerr<<"Two points coincide.\n";
					return false;
				}
			}
		for(int i=0;i<n;i++)
			for(int j=0;j<n;j++)
				for(int k=0;k<n;k++){
					set<int> pts={i,j,k};
					if(pts.size()!=3) continue;
					if(on_the_line(pt[i],pt[j],pt[k])){
						if (veb)std::cerr<<"Three points which are colinear.\n";
						return false;
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
				if (veb)std::cerr<<"n="<<n<<"\n";
				if (veb)for(int k=0;k<n;k++) std::cerr<<"   Point #"<<k+1<<": ("<<pt[k].x<<","<<pt[k].y<<")\n";
				if (veb)std::cerr<<"A 6-hole:\n";
				if (veb)for(int i=0;i<6;i++) std::cerr<<"   "<<a[i]<<": ("<<pt[a[i]].x<<","<<pt[a[i]].y<<")\n";
				return false;
			}
		}
		return true;
	}
}

vector<Point> shake_method(vector<Point> pt,LL T){
	static const LL dx[4] = {-1, 0, 0, 1};
	static const LL dy[4] = {0, -1, 1, 0};
	vector<Point> tmp;
	tmp = pt;
	for(int i=0;i<pt.size()-1;i++){
		int chosen = randint(-4, 3);
		if(chosen>=0){
			tmp[i].x += (1<<T)*dx[chosen];
			tmp[i].y += (1<<T)*dy[chosen];
		}
	}
	return tmp;
} 

vector<Point> shake(vector<Point> pt, LL T = 0) {
	if(T<0) return pt;
	vector<Point> tmp;
	for(;T>=0;T--){
		tmp = shake_method(pt,T);
		if (ahdoc::check_no6hole(tmp))
			pt=tmp;
	}
	return pt;
}

int main(){
	vector<vector<Point>> in={{
{11713,13448},{-191712,-38364},{-8504,22138},{4619,-35177},{92495,-39193},{-133,12964},{9047,2634},{204756,-195599},{-7099,8151},{-18223,-45783},{5238,101447},{14701,141891},{8914,108611},{-122945,-19363},{-398,-1207},{372,124862},{-32489,66505},{6059,25295},{-2465,28887},{-84239,39230},{-27893,-20741},{-69969,27293},{43326,29383},{26397,226741},{73332,-33331},{120134,-93838}
},{
{25804,-42493},{147227,-202204},{1804,141},{12476,-13424},{-20789,217372},{4027,54369},{-32178,100498},{36392,22009},{27347,-25192},{-969,-27059},{-10020,5381},{13857,2190},{-17610,133067},{19540,19887},{-4516,-10046},{-123788,17888},{20008,-6181},{79559,-39129},{-236216,43261},{52620,-14994},{-4788,37188},{179970,-79153},{354,33475},{-87438,-47007},{-207946,-181231},{-102449,5049}
},{
{-1317,625},{-57290,41809},{5770,-46080},{27957,9921},{207267,-175700},{-10251,-4343},{22440,34550},{15871,213},{132987,145979},{-24889,53226},{-16996,20990},{7712,-7515},{59539,3709},{72452,78037},{-23349,10828},{77,5514},{-218587,-220829},{-105608,38676},{-129106,162371},{24881,-10225},{-232366,200029},{-38382,8525},{18893,16312},{63464,-33301},{-96180,-96300},{50406,-61676}
},{
{10029,47390},{102191,-35525},{-23547,-63602},{-20069,-9386},{230955,-85771},{-303,-2262},{26952,16962},{-223716,205391},{-9433,-14039},{-26045,17391},{35891,11985},{28020,204664},{55500,55038},{-93551,22211},{-3846,8488},{-119793,-128815},{1926,-22401},{-89600,-64488},{56018,-12495},{68943,6643},{18262,-1276},{21114,-10650},{-24724,-25869},{-204214,-210578},{76720,73408},{242557,243898}
},{
{-16180,-2460},{70432,120666},{7127,-9744},{-1505,282},{58859,-97836},{-37881,47880},{94631,-213428},{23069,-9673},{58701,-125389},{132835,226180},{-15815,62465},{-86856,222899},{68465,-48141},{24489,-68770},{1611,7792},{-10986,19084},{19512,26940},{-16162,157789},{-9989,-46916},{-7142,74207},{-3381,16357},{-34100,86966},{2162,42995},{-208384,-144332},{-100591,89706},{-223141,195629}
},{
{-132600,-26428},{22913,-111245},{8752,-30276},{223,12679},{-6256,-2194},{-16423,-15083},{-1955,2728},{-41776,23540},{-175248,-7076},{7385,-35374},{145408,-227520},{-103792,-53303},{15857,-17637},{-218040,1498},{8381,21084},{-8004,12561},{22834,-5566},{-154046,-14600},{-19952,2644},{66249,-147944},{22853,-43124},{-22739,66678},{-31209,122339},{-8584,43292},{44540,196951},{18583,95147}
},{
{5823,-8117},{-28119,7045},{8468,-108048},{-4442,40751},{-96639,-98717},{-6,-10239},{53753,-29228},{-1076,267},{741,-46181},{-8229,-53059},{8272,-156082},{-15306,-119498},{19680,-231579},{1805,-13198},{9224,-66472},{-8399,-18880},{9283,63933},{23131,-106985},{78041,202391},{2194,3160},{27482,-40438},{9263,17843},{227019,-100953},{-39661,157550},{-192993,-232984},{-55388,227690}
}};
	
	int T=1000;
	for(auto pt:in){
		int cnt=0;
		for(int i=0;i<T;i++){
			vector<Point> tmp;
			tmp = shake_method(pt,7);
			if (ahdoc::check_no6hole(tmp))
				++cnt;
		}
		cerr<<cnt<<" out of "<<T<<" successful shakings!\n";
	}
} 


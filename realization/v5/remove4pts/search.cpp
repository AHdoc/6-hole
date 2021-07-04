#include<iostream>
#include<fstream>
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

#include <chrono>
using namespace std::chrono;
long long gettime() {
	return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

long long randint(long long l, long long r) {
	static mt19937_64 generator(123);
	return uniform_int_distribution<long long>(l, r)(generator);
	// return rand() % (r - l + 1) + l;
}

typedef long long LL;

namespace ahdoc{
	typedef long long LL;
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
long long xxtime;

void print_log() {
	cerr<<"time:"<<(gettime() - xxtime)<<endl;
	cerr<<"tot_6holes:"<<tot_b<<endl;
	cerr<<"tot_query_find6hole = "<<tot_query_find6hole<<" & tot_query_check = "<<tot_query_check<<":";
	for(int j=1;j<=30;j++){
		if(achievement[j]==0) continue;
		cerr<<j<<":"<<achievement[j]<<" ";
	}
	cerr<<"\n";

}

bool find6hole(vector<pair<LL,LL>> pt,pair<LL,LL> p){
	bool ret=ahdoc::find6hole(pt,p);
	++tot_query_find6hole;
	tot_b+=ret;
	if(tot_query_find6hole%1000000LL==0){
		print_log();
	}
	return ret;
}


bool check(const vector<pair<LL,LL>> &pt, const pair<LL,LL> &p){
	++tot_query_check;
	set<pair<LL,LL>> ang;
	for(int j=1;j<=pt.size();j++){
		LL x=pt[j-1].x-p.x,y=pt[j-1].y-p.y;
		if(x==0 || y==0) return false;
		LL g=gcd(absLL(x),absLL(y)); x/=g; y/=g;
		if(x<0){x=-x; y=-y;}
		else if(x==0 && y<0) y=-y;
		ang.insert({x,y}); 
	}
	if(ang.size()!=pt.size()) return false;
	return !find6hole(pt,p);
}


bool check_no6hole(vector<pair<LL,LL>> pt){ // Check whether or not the given set of points contains no 6-hole
	vector<pair<LL,LL>> pt2;
	pt2.clear();
	pt2.push_back(pt[0]);
	for(int i=1;i<pt.size();i++)
		if(check(pt2,pt[i]))
			pt2.push_back(pt[i]);
		else
			return false;
	return true;
}

bool check_no6hole_bruteforce(vector<pair<LL,LL>> pt, bool veb=false){ // Check whether or not the given set of points contains no 6-hole
	int n=pt.size();
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++){
			set<int> pts={i,j};
			if(pts.size()!=2) continue;
			if(pt[i].x==pt[j].x && pt[i].y==pt[j].y){
				if (veb)cerr<<"Two points coincide.\n";
				return false;
			}
		}
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			for(int k=0;k<n;k++){
				set<int> pts={i,j,k};
				if(pts.size()!=3) continue;
				if(on_the_line(pt[i],pt[j],pt[k])){
					if (veb)cerr<<"Three points which are colinear.\n";
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
			if (veb)cerr<<"n="<<n<<"\n";
			if (veb)for(int k=0;k<n;k++) cerr<<"   Point #"<<k+1<<": ("<<pt[k].x<<","<<pt[k].y<<")\n";
			if (veb)cerr<<"A 6-hole:\n";
			if (veb)for(int i=0;i<6;i++) cerr<<"   "<<a[i]<<": ("<<pt[a[i]].x<<","<<pt[a[i]].y<<")\n";
			return false;
		}
	}
	return true;
}

LL totdfs=0;
LL N = 1000000, M = 16*N, T=100;
ofstream fout("tmp_a.txt");
LL largest = 0;

void add_pts(const vector<pair<LL, LL>> &pt) {
	if (achievement[pt.size()] == 0) {
		largest = pt.size();
	}
	if (pt.size() >= largest && achievement[pt.size()] <= 100) {
		fout << pt.size() << ":" << endl;
		for(int i=1;i<=pt.size();i++) 
			fout<<"("<<pt[i-1].x<<","<<pt[i-1].y<<")";
		fout<<endl;
	}
	achievement[pt.size()]++;
	
	if(pt.size()>=30){
		for(int i=1;i<=pt.size();i++) cerr<<"Point #"<<i<<": ("<<pt[i-1].x<<","<<pt[i-1].y<<")\n";
		if (!check_no6hole_bruteforce(pt, true))
			exit(1);
		exit(0);
	}
}

const LL dx[4]={-1,0,0,1};
const LL dy[4]={0,-1,1,0};

template<typename T>
vector<T> sample(const vector<T> &vec, int n) {
	vector<int> ids(vec.size());
	for (int i = 0; i < vec.size(); i++)
		ids[i] = i;
	random_shuffle(ids.begin(), ids.end());
	vector<T> ret(n);
	for (int i = 0; i < n; i++)
		ret[i] = vec[ids[i]];
	return ret;
}

vector<pair<LL, LL>> remove_et_shake(vector<pair<LL, LL>> pt,int steps=0) {
	int n=pt.size();
	vector<pair<LL, LL>> tmp;
	
	for(int u,v;;){ 
		tmp = sample(pt, pt.size() - 3);
		if(check_no6hole(tmp)){
			pt=tmp;
			break;
		}
	}
	
	/*
	n=pt.size();
	tmp.resize(n);
	LL R = 1, iter;
	for (int step = 0; step < steps; step++) {
		for (iter = 0; iter < 1000; iter++) {
			for (int i = 0; i < n; i++) {
				tmp[i] = pt[i];
				int chosen=randint(0, 3);
				tmp[i].x += dx[chosen];
				tmp[i].y += dy[chosen];
			}
			if (check_no6hole(tmp)) {
				pt = tmp;
				break;
			}
		}
	}*/
	return pt;
}

vector<pair<LL, LL>> remove(vector<pair<LL, LL>> pt, int D = 3) {
	vector<pair<LL, LL>> tmp;
	while (1) {
		tmp = sample(pt, pt.size() - D);
		if (check_no6hole(tmp))
			return tmp;
	}
	return vector<pair<LL, LL>>();
}

vector<pair<LL, LL>> shake(vector<pair<LL, LL>> pt, LL R = 1, LL T = 1000) {
	static const LL dx[4] = {-1, 0, 0, 1};
	static const LL dy[4] = {0, -1, 1, 0};
	vector<pair<LL, LL>> tmp(pt.size());
	for (int iter = 0; iter < T; iter++) {
		for (int i = 0; i < n; i++) {
			tmp[i] = pt[i];
			int chosen = randint(0, 3);
			tmp[i].x += dx[chosen];
			tmp[i].y += dy[chosen];
		}
		if (check_no6hole(tmp))
			return tmp;
	}
	return vector<pair<LL, LL>>();
}
vector<pair<LL, LL>> expand(vector<pair<LL, LL>> pt, LL R = 5000, LL T = 1000) {
	pair<LL,LL> p;
	while (T--) {
		p.x = randint(-R, R);
		p.y = randint(-R, R);
		if (check(pt, p)) {
			pt.push_back(p);
			return pt;
		}
	}
	return vector<pair<LL, LL>>();
}

void dfs(vector<pair<LL,LL>> pt) {
	cerr << "dfs:" << endl;
	vector<pair<LL,LL>> tmp;
	while(1) {
		add_pts(pt);
		totdfs++;
		
		LL R = 500 * pt.size();
		LL T = 1000 + max((LL)pt.size() - 20, 1LL)*max((LL)pt.size() - 20, 1LL) * 1000;// + max((LL)pt.size() - 22, 0LL) * 10000;

		tmp = expand(pt, R, T);
		if (tmp.empty()) {
			if (pt.size() >= 22)
				pt=remove_et_shake(pt , 100);
			else{
				pt.clear();
				pt.push_back({ 0, 0 });
			}
		}
		else
			pt = tmp;
	}
}

struct BeamSearch {
	struct State {
		vector<pair<LL, LL>> pt;
	};
	vector<queue<State>> Q;
	LL head, width;

	BeamSearch(int limit = 31, LL width = 1000) : width(width) {
		Q.resize(limit);
		head = 0;
		vector<pair<LL, LL>> tmp = {{0, 0}};
		Add(tmp);
	}

	void log() {
		static LL cnt = 0;
		if (++cnt % 10000 == 0 && head >= 1) {
			cerr << "Q: ";
			for (int i = 1; i <= head; i++)
				if (!Q[i].empty()) 
					cerr << i<<":"<<Q[i].size() << " ";
			cerr << endl;
		}
	}

	void Add(const vector<pair<LL, LL>> &pt) {
		if (!pt.empty()) {
			if (Q[pt.size()].size() < width)
				Q[pt.size()].push(State({ pt }));
			head = max(head, (LL)pt.size());
			add_pts(pt);
			log();
		}
	}

	State pop() {
		while (head > 0 && Q[head].empty())
			head--;
		State ret = Q[head].front();
		if (head > 1)
			Q[head].pop();
		return ret;
	}

	void run() {
		while (1) {
			const State &stat = pop();
			const vector<pair<LL, LL>> &pt = stat.pt;
			// add
			LL R = 10000 * pt.size();
			LL T = 1000 + max((LL)pt.size() - 20, 0LL)*max((LL)pt.size() - 20, 0LL) * 1000;
			LL add_iters = 1;
			for (int iter = 0; iter < add_iters; iter++) {
				const vector<pair<LL, LL>> &tmp = expand(pt, R, T);
				Add(tmp);
				// if (tmp.empty() && pt.size() >= 22) {
				// 	const vector<pair<LL, LL>> &tmp = remove(pt);
				// 	Add(tmp);
				// }
			}
			// del
			LL del_iters = 0;
			if (pt.size() >= 21) {
				del_iters = 2;
				// del_iters = (pt.size() - 20);
				// del_iters *= del_iters * del_iters;
			}
			for (int iter = 0; iter < del_iters; iter++) {
				const vector<pair<LL, LL>> &tmp = remove(pt, min(4, iter + 3));
				Add(tmp);
			}
			//shake
			// LL shake_iters = 0;
			// if (pt.size() > 21) {
			// 	shake_iters = (pt.size() - 20);
			// 	shake_iters *= shake_iters * shake_iters;
			// }
			// for (int iter = 0; iter < shake_iters; iter++) {
			// 	const vector<pair<LL, LL>> &tmp = shake(pt, 1);
			// 	Add(tmp);
			// }
		}
	}
};

void Realizer(){
	// vector<pair<LL,LL>> pt;
	// pt.push_back({ 0, 0 });
	// for(;;) dfs(pt); // must fix the kernel point
	BeamSearch bm(31, 100);
	bm.run();
}

int main() {
	// srand(time(NULL));
	srand(123);
	xxtime = gettime();
	Realizer();
}

#include <iostream>
#include <fstream>
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
#define x first
#define y second

bool atan2_less(const pair<LL,LL> &p1, const pair<LL,LL> &p2) {
	// return atan2(p1.y,p1.x)<atan2(p2.y,p2.x);
	int s1 = (((p1.y > 0) << 1) - 1) << (p1.x < 0);
	int s2 = (((p2.y > 0) << 1) - 1) << (p2.x < 0);
	if (s1 != s2)
		return s1 < s2;
	return (p1.x * p2.y - p2.x * p1.y) > 0;
}

namespace ahdoc{
	typedef long long LL;
	#define x first
	#define y second
	#define c1 first.first
	#define c2 first.second
	#define c3 second
	#define MP3(a,b,c) make_pair(make_pair(a,b),c)
	
	const int MAXN=30;
	
	thread_local int n;
	thread_local pair<LL,LL> pt[MAXN];
	thread_local int Q_head[MAXN],Q_tail[MAXN];
	thread_local pair<pair<int,int>,int> Q[MAXN][MAXN];
	thread_local pair<pair<int,int>,int> C[MAXN];
	
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
			//	cerr<<"	   Q["<<k<<"]={ ";
			//	for(auto x:Q[k]) cerr<<"("<<x.c1<<","<<x.c2<<","<<x.c3<<") ";
			//	cerr<<"}\n";
			//	cerr<<"	   C["<<k<<"]=("<<C[k].c1<<","<<C[k].c2<<","<<C[k].c3<<")\n";
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
		//sort(pt,pt+n,[](pair<LL,LL> p1,pair<LL,LL> p2){return atan2(p1.y,p1.x)<atan2(p2.y,p2.x);});
		sort(pt,pt+n,atan2_less);
		if(solve()) return true;
		for(int i=0;i<n;i++){
			pt[i].x=p.x-_pt[i].x;
			pt[i].y=p.y-_pt[i].y;
		}
		//sort(pt,pt+n,[](pair<LL,LL> p1,pair<LL,LL> p2){return atan2(p1.y,p1.x)<atan2(p2.y,p2.x);});
		sort(pt,pt+n,atan2_less);
		if(solve()) return true;
		return false;
	}
}



LL absLL(LL x){return (x<0?-x:x);}
LL gcd(LL a,LL b){return (b==0?a:gcd(b,a%b));}

const int MAXN=30;

int n;
LL tot_query_check,tot_query_find6hole,achievement[MAXN+1];
LL radius[MAXN+1],lvl[MAXN+2];

LL crossproduct(const pair<LL,LL> &a,const pair<LL,LL> &b,const pair<LL,LL> &c){return (b.x-a.x)*(c.y-a.y)-(c.x-a.x)*(b.y-a.y);}
LL innerproduct(const pair<LL,LL> &a,const pair<LL,LL> &b,const pair<LL,LL> &c){return (b.x-a.x)*(c.x-a.x)+(b.y-a.y)*(c.y-a.y);}
LL crossproduct(const pair<LL,LL> &a,const pair<LL,LL> &b,const pair<LL,LL> &c,const pair<LL,LL> &d){return (b.x-a.x)*(d.y-c.y)-(d.x-c.x)*(b.y-a.y);}
bool on_the_left(const pair<LL,LL> &a,const pair<LL,LL> &b,const pair<LL,LL> &c){return crossproduct(a,b,c)>0;}
bool on_the_line(const pair<LL,LL> &a,const pair<LL,LL> &b,const pair<LL,LL> &c){return crossproduct(a,b,c)==0;}

long long tot_b=0; 
long long xxtime;

string xxpslog();

void print_log() {
	cerr<<"time:"<<(gettime() - xxtime)<<endl;
	cerr<<"tot_6holes:"<<tot_b<<endl;
	cerr<<"tot_query_find6hole = "<<tot_query_find6hole<<" & tot_query_check = "<<tot_query_check<<":";
	for(int j=1;j<=30;j++){
		if(achievement[j]==0) continue;
		cerr<<j<<":"<<achievement[j]<<" ";
	}
	cerr << " @ " << xxpslog();
	cerr<<"\n";

}

bool find6hole(vector<pair<LL,LL>> pt,pair<LL,LL> p){
	bool ret=ahdoc::find6hole(pt,p);
	#pragma omp critical
	{
		++tot_query_find6hole;
		tot_b+=ret;
		if(tot_query_find6hole%10000000LL==0){
			print_log();
		}
	}
	return ret;
}


bool check(const vector<pair<LL,LL>> &pt, const pair<LL,LL> &p){
	#pragma omp critical
	++tot_query_check;
	// set<pair<LL,LL>> ang;
	// for(int j=1;j<=pt.size();j++){
	// 	LL x=pt[j-1].x-p.x,y=pt[j-1].y-p.y;
	// 	if(x==0 || y==0) return false;
	// 	LL g=gcd(absLL(x),absLL(y)); x/=g; y/=g;
	// 	if(x<0){x=-x; y=-y;}
	// 	else if(x==0 && y<0) y=-y;
	// 	ang.insert({x,y}); 
	// }
	// if(ang.size()!=pt.size()) return false;
	for (int i = 0; i < pt.size(); i++) {
		LL x = pt[i].x - p.x, y = pt[i].y - p.y;
		if (x == 0 || y == 0)
			return false;
		for (int j = 0; j < i; j++)
			if (x * (pt[j].y - p.y) == y * (pt[j].x - p.x))
				return false;
	}
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


class Serializer {
public:
	std::fstream file;
public:
	Serializer() {}
	Serializer(const std::string &filename) : file(filename, std::ios::out | std::ofstream::binary) {}

	void close() {
		file.close();
	}

	Serializer& operator <<(long long t) {
		file.write((char *)&t, sizeof(t));
		return *this;
	}

	template<typename T>
	Serializer& operator <<(std::vector<T> &vec) {
		(*this) << (long long)vec.size();
		for (T &t : vec)
			(*this) << t;
		return *this;
	}

	template<typename T1, typename T2>
	Serializer& operator <<(std::pair<T1, T2> &t) {
		return (*this) << t.first << t.second;
	}

	template<typename T>
	Serializer& operator <<(T &t) {
		return t.Serialize(*this);
	}
};

class Deserializer {
public:
	std::fstream file;
public:
	Deserializer() {}
	Deserializer(const std::string &filename) : file(filename, std::ios::in | std::ofstream::binary) {}

	void close() {
		file.close();
	}

	Deserializer& operator >>(long long &t) {
		file.read((char *)&t, sizeof(t));
		return *this;
	}

	template<typename T>
	Deserializer& operator >>(std::vector<T> &vec) {
		long long n;
		(*this) >> n;
		vec.resize(n);
		for (size_t i = 0; i < n; i++)
			(*this) >> vec[i];
		return *this;
	}

	template<typename T1, typename T2>
	Deserializer& operator >>(std::pair<T1, T2> &t) {
		return (*this) >> t.first >> t.second;
	}

	template<typename T>
	Deserializer& operator >>(T &t) {
		return t.Deserialize(*this);
	}
};

LL totdfs=0;
LL N = 1000000, M = 16*N, T=100;
ofstream fout;
LL largest = 0;

void log_pts(const vector<pair<LL, LL>> &pt) {
	fout << pt.size() << ":" << endl;
	for(int i=1;i<=pt.size();i++) 
		fout<<"("<<pt[i-1].x<<","<<pt[i-1].y<<")";
	fout<<endl;
}

void add_pts(const vector<pair<LL, LL>> &pt) {
	if (achievement[pt.size()] == 0) {
		largest = pt.size();
	}
	// if (pt.size() >= largest && achievement[pt.size()] <= 100) 
	// 	log_pts(pt);
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

/*vector<pair<LL, LL>> remove_et_shake(vector<pair<LL, LL>> pt,int steps=0) {
	int n=pt.size();
	vector<pair<LL, LL>> tmp;
	
	for(int u,v;;){ 
		tmp = sample(pt, pt.size() - 3);
		if(check_no6hole(tmp)){
			pt=tmp;
			break;
		}
	}
	
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
	}
	return pt;
}*/

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
			pt = tmp;
			//return tmp;
	}
	return pt;
	//return vector<pair<LL, LL>>();
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
/*
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
}*/

struct BeamSearch {
	struct State {
		vector<pair<LL, LL>> pt;
		// shared_ptr<State> pre;
		// std::string m;
		// State() {}
		// State(const vector<pair<LL, LL>> &pt = vector<pair<LL, LL>(), shared_ptr<State> pre = nullptr, const std::string &m) : pt(pt), pre(pre), m(m) {}
	};
	Serializer& Serialize(Serializer &ar) {
		ar << head << width << largest << largest_pt;
		for (LL i = 0; i < Q.size(); i++) {
			vector<vector<pair<LL, LL>>> tmp;
			while (!Q[i].empty()) {
				tmp.push_back(Q[i].front().pt);
				Q[i].pop();
			}
			ar << tmp;
			for (LL j = 0; j < tmp.size(); j++)
				Q[i].push({tmp[j]});
		}
		return ar;
	}
	Deserializer& Deserialize(Deserializer &ar) {
		ar >> head >> width >> largest >> largest_pt;
		for (LL i = 0; i < Q.size(); i++) {
			vector<vector<pair<LL, LL>>> tmp;
			ar >> tmp;
			for (LL j = 0; j < tmp.size(); j++)
				Q[i].push({tmp[j]});
		}
		return ar;
	}
	vector<queue<State>> Q;
	LL head, width, largest;
	vector<pair<LL, LL>> largest_pt;

	BeamSearch(int limit = 31, LL width = 1000) : width(width) {
		Q.resize(limit);
		head = 0;
		largest = 0;
		vector<pair<LL, LL>> tmp = {{0, 0}};
		// Add(State({tmp, nullptr, ""}));
		Add(State({tmp}));
	}

	// void log() {
	// 	thread_local LL cnt = 0;
	// 	if (++cnt % 10000 == 0 && head >= 1) {
	// 		cerr << "Q: ";
	// 		for (int i = 1; i <= head; i++)
	// 			if (!Q[i].empty()) 
	// 				cerr << i<<":"<<Q[i].size() << " ";
	// 		cerr << endl;
	// 	}
	// }
	
	// void log_all(State stat) {
	// 	log_pts(stat.pt);
	// 	fout << "by " << stat.m << " from" << endl;
	// 	shared_ptr<State> now = stat.pre;
	// 	while (1) {
	// 		log_pts(now->pt);
	// 		if (now->pt.size() <= 1)
	// 			break;
	// 		fout << "by " << now->m << " from" << endl;
	// 		now = now->pre;
	// 	}
	// }

	void Add(const State &stat) {
		if (!stat.pt.empty()) {
			if (Q[stat.pt.size()].size() < width)
				Q[stat.pt.size()].push(stat);
			if (largest < stat.pt.size()) {
				// log_pts(pt);
				// if (stat.pt.size() >= 25)
					// #pragma omp critical
					// log_all(stat);
				largest = stat.pt.size();
				largest_pt = stat.pt;
			}
			head = max(head, (LL)stat.pt.size());
			#pragma omp critical
			{
				add_pts(stat.pt);
				//#pragma omp critical
				// log();
			}
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

	void step() {
		const State &stat = pop();
		// shared_ptr<State> pre = make_shared<State>(stat);
		const vector<pair<LL, LL>> &pt = stat.pt;
		// add
		LL R = 100 * pt.size();
		// double Rr = 10;
		// for (int q = pt.size() - 1; q > 0; q--) Rr *= 1.26;
		// LL R = (LL)(Rr + 1);
		// if (pt.size() == largest) cerr << pt.size() << " : " << R << endl; 
		// LL R = 1000;
		LL T = (pt.size()<=20) ? (3LL*(LL)pt.size()*(LL)pt.size()) : (1000 + max((LL)pt.size() - 20, 0LL)*max((LL)pt.size() - 20, 0LL) * 1000);
		LL add_iters = 1;
		for (int iter = 0; iter < add_iters; iter++) {
			const vector<pair<LL, LL>> &tmp = expand(pt, R, T);
			// Add(State({tmp, pre, "add"}));
			Add(State({tmp}));
			// if (tmp.empty() && pt.size() >= 22) {
			// 	const vector<pair<LL, LL>> &tmp = remove(pt);
			// 	Add(tmp, pre);
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
			int c = 3;//min(3, iter + 2);
			const vector<pair<LL, LL>> &tmp = shake(remove(pt, c),1,100);
			// Add(tmp, pre);
			// Add(State({tmp, pre, "del_" + std::to_string(c)}));
			Add(State({tmp}));
		}
		//shake
		// LL shake_iters = 0;
		// if (pt.size() > 21) {
		// 	shake_iters = (pt.size() - 20);
		// 	shake_iters *= shake_iters * shake_iters;
		// }
		// for (int iter = 0; iter < shake_iters; iter++) {
		// 	const vector<pair<LL, LL>> &tmp = shake(pt, 1);
		// 	Add(tmp, pre);
		// }
	}
	void run() {
		while (1)
			step();
	}
};

/*
struct Pool {
	struct State {
		vector<pair<LL, LL>> pt;
		shared_ptr<State> pre;
		std::string m;
		// State() {}
		// State(const vector<pair<LL, LL>> &pt = vector<pair<LL, LL>(), shared_ptr<State> pre = nullptr, const std::string &m) : pt(pt), pre(pre), m(m) {}
	};
	vector<State> pool;
	LL head, width, largest;

	Pool(LL width = 1000) : width(width) {
		head = 0;
		largest = 0;
		vector<pair<LL, LL>> tmp = {{0, 0}};
		Add(State({tmp, nullptr, ""}));
	}
	
	void log_all(State stat) {
		log_pts(stat.pt);
		fout << "by " << stat.m << " from" << endl;
		shared_ptr<State> now = stat.pre;
		while (1) {
			log_pts(now->pt);
			if (now->pt.size() <= 1)
				break;
			fout << "by " << now->m << " from" << endl;
			now = now->pre;
		}
	}

	void Add(const State &stat) {
		if (!stat.pt.empty()) {
			if (largest < stat.pt.size()) {
				// log_pts(pt);
				if (stat.pt.size() >= 25) {
					#pragma omp critical
					log_all(stat);
				}
				largest = stat.pt.size();
				pool.clear();
				pool.push_back(stat);
				#pragma omp critical
				add_pts(stat.pt);
			}
			else if (largest == stat.pt.size()) {
				if (pool.size() < width)
					pool.push_back(stat);
				else
					pool[randint(0, pool.size() - 1)] = stat;
				#pragma omp critical
				add_pts(stat.pt);
			}
		}
	}

	State pop() {
		return pool[randint(0, pool.size() - 1)];
	}

	vector<pair<LL, LL>> expands_xx(vector<pair<LL, LL>> pt) {
		vector<pair<LL, LL>> tmp;
		LL init = pt.size();
		while (1) {
			LL n = pt.size();
			LL R = 100 * n, T = 1000 + max((LL)n - 20, 0LL)*max((LL)n - 20, 0LL) * 1000;
			tmp = expand(pt, R, T);
			if (tmp.empty())
				break;
			pt = tmp;
		}
		if (pt.size() <= init)
			return vector<pair<LL, LL>>();
		return pt;
	}

	void step() {
		const State &stat = pop();
		shared_ptr<State> pre = make_shared<State>(stat);
		const vector<pair<LL, LL>> &pt = stat.pt;
		// add
		LL add_iters = 1;
		for (int iter = 0; iter < add_iters; iter++) {
			const vector<pair<LL, LL>> &tmp = expands_xx(pt);
			Add(State({tmp, pre, "add"}));
		}
		// del
		LL del_iters = 0;
		if (pt.size() >= 21) {
			del_iters = 2;
		}
		for (int iter = 0; iter < del_iters; iter++) {
			int c = min((int)pt.size() - 1, iter + 2);
			vector<pair<LL, LL>> tmp = remove(pt, c);
			tmp = expands_xx(tmp);
			Add(State({tmp, pre, "del_" + std::to_string(c)}));
		}
		//shake
		// LL shake_iters = 0;
		// if (pt.size() > 21) {
		// 	shake_iters = (pt.size() - 20);
		// 	shake_iters *= shake_iters * shake_iters;
		// }
		// for (int iter = 0; iter < shake_iters; iter++) {
		// 	const vector<pair<LL, LL>> &tmp = shake(pt, 1);
		// 	Add(tmp, pre);
		// }
	}
	void run() {
		while (1)
			step();
	}
};*/

LL nn = 1000, w = 10000;
vector<BeamSearch> ps(nn, BeamSearch(31, w));
// vector<Pool> ps(nn, Pool(w));

string xxpslog() {
	const int N = 40;
	LL cnt[N];
	for (int i = 0; i < N; i++)
		cnt[i] = 0;
	for (int i = 0; i < nn; i++)
		cnt[ps[i].largest]++;
	string ret = "";
	for (int i = 0; i < N; i++)
		if (cnt[i] > 0) {
			if (ret != "")
				ret += " ";
			ret += to_string(i) + ":" + to_string(cnt[i]);
		}
	return ret;
}

std::string log_filename = "tmp_a.txt";
std::string cache_filename = "tmp";
std::string start_filename = "";//"tmp1.binary";
void Realizer(){
	// vector<pair<LL,LL>> pt;
	// pt.push_back({ 0, 0 });
	// for(;;) dfs(pt); // must fix the kernel point

	if (start_filename != "") {
		cerr << "start loading from " << start_filename << endl;
		Deserializer(start_filename) >> ps;
		cerr << "progress loaded from " << start_filename << " @ " << xxpslog() << endl;
	}
	
	LL log_time = gettime(), cache_id = 0, interval = 1 * 60 * 1000;
	while (1) {
		#pragma omp parallel for
		for (int i = 0; i < ps.size(); i++)
			ps[i].step();
		if ((gettime() - log_time) >= interval) { //1h
			string filename = cache_filename + "_" + to_string(cache_id % 2 + 1) + ".binary";
			cerr << "start saving to " << filename << " @ " << xxpslog() << endl;
			Serializer(filename) << ps;
			cerr << "progress saved to " << filename << endl;
			cache_id++;
			log_time = gettime();
		}
	}
}

int main(int argc, char *argv[]) {
	// std::ios::sync_with_stdio(false);
	if (argc > 1) {
		auto peek = [&](int idx) {
			if (idx < argc)
				return std::string(argv[idx]);
			else {
				std::cerr << "error in augment " << idx << std::endl;
				exit(2);
			}
		};
		for (int i = 1; i < argc; i++)
			if (argv[i][0] == '-') {
				char cmd = argv[i][1];
				if (cmd == 'o')
					log_filename = peek(i + 1);
				else if (cmd == 'c')
					cache_filename = peek(i + 1);
				else if (cmd == 's')
					start_filename = peek(i + 1);
			}
	}
	fout.open(log_filename);
	// srand(time(NULL));
	srand(123);
	omp_set_num_threads(7);
	xxtime = gettime();
	Realizer();
}

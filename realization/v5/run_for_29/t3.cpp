/*
New shaking method.
*/

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
#include <functional>
#include <algorithm>
#include <unordered_map>
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

	bool atan2_less(Point p1, Point p2) {
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

	bool check_linear(Point pt[], int n) {
		// given pt sorted in polar order, check whether no two points being colinear with the origin
		for(int i=1;i<n;i++)
			if(pt[i-1].x*pt[i].y==pt[i].x*pt[i-1].y)
				return false;
		for(int i=0,j=1;i<n;i++){
			while(pt[i].x*pt[j].y-pt[j].x*pt[i].y>0) j=(j+1)%n;
			if(i!=j && pt[i].x*pt[j].y==pt[j].x*pt[i].y)
				return false;
		}
		return true;
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
		if (n >= 2 && !check_linear(pt, n))
			return false;
		// if(n>=2){
		// 	for(int i=1;i<n;i++)
		// 		if(pt[i-1].x*pt[i].y==pt[i].x*pt[i-1].y)
		// 			return false;
		// 	for(int i=0,j=1;i<n;i++){
		// 		while(pt[i].x*pt[j].y-pt[j].x*pt[i].y>0) j=(j+1)%n;
		// 		if(i!=j && pt[i].x*pt[j].y==pt[j].x*pt[i].y)
		// 			return false;
		// 	}
		// }

		
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

string log_beams();

struct Config {
	LL n, w, t, b, k, r, lt;
	LL cache_interval, log_interval;
	string log_filename, cache_filename, start_cache_filename, start_log_filename;
	vector<int> dels;

	std::unordered_map<std::string, std::string> c, args;
	std::vector<std::string> helps;
	std::vector<std::function<void()>> fns;

	void add(const std::string &id, const std::string &t, const std::string &help, const std::string &dv) {
		std::string hs = "    -" + id;
		if (!t.empty())
			hs += " " + t;
		hs += ": " + help;
		if (!dv.empty())
			hs += "[default:" + dv + "]";
		helps.push_back(hs);
		if (!c.count(id))
			c[id] = dv;
		if (!t.empty())
			args[id] = t;
	}
	template<typename T>
	void add_int(T &v, const std::string &id, const std::string &t, const std::string &help, const std::string &dv) {
		add(id, t, help, dv);
		fns.push_back([&v, id, this](){
			v = std::stoll(c[id]);
		});
	}
	template<typename T>
	void add_str(T &v, const std::string &id, const std::string &t, const std::string &help, const std::string &dv) {
		add(id, t, help, dv);
		fns.push_back([&v, id, this](){
			v = c[id];
		});
	}
	template<typename T>
	void add_vec(T &v, const std::string &id, const std::string &t, const std::string &help, const std::string &dv) {
		add(id, t, help, dv);
		fns.push_back([&v, id, this](){
			v = T(parse_vec(c[id]));
		});
	}

	Config() {
		add_str(log_filename, "o", "string", "output file", "");
		add_str(cache_filename, "c", "string", "cache file, no cache if empty", "");
		add_str(start_log_filename, "l", "string", "init log file, no init loading if empty", "");
		add_str(start_cache_filename, "s", "string", "init cache file, no init loading if empty", "");
		add_int(n, "n", "integer", "number of beams", "3000");
		add_int(b, "b", "integer", "number of least loaded pts", "1");
		add_int(k, "k", "integer", "number of shakes", "10");
		add_int(t, "t", "integer", "number of threads", "7");
		add_int(r, "r", "integer", "size of radius coefficient", "10000");
		add_int(w, "w", "integer", "size of beams", "1000");
		add_int(lt, "lt", "integer", "keep only points of this size from init log file", "-1");
		add_vec(dels, "d", "list", "number points to delete", "3,3");
		add_int(cache_interval, "ci", "integer", "number of seconds between cache", "3600");
		add_int(log_interval, "li", "integer", "max number of seconds between logs", "300");
		add("h", "", "print help message and quit", "");
	}

	string to_string() {
		string ret;
		for (auto p : args) {
			if (ret != "")
				ret += ", ";
			ret += p.first + "=" + c[p.first];
		}
		return ret;
	}

	vector<int> parse_vec(const std::string &s) {
		vector<int> ret;
		string tmp;
		for (int i = 0, j = 0; i <= s.length(); i++) {
			if (i == s.length() || s[i] == ',') {
				ret.push_back(stoi(tmp));
				tmp.clear();
			} else
				tmp += s[i];
		}
		return ret;
	}

	void parse(int argc, char *argv[]) {
		auto peek = [&](int idx) {
			if (idx < argc)
				return std::string(argv[idx]);
			else {
				std::cerr << "error in augment " << idx << std::endl;
				print_help();
				exit(2);
			}
		};
		for (int i = 1; i < argc; i++)
			if (argv[i][0] == '-') {
				string cmd = string(argv[i] + 1);
				if (cmd == "h")
					print_help();
				else {
					std::string s = peek(i + 1);
					c[cmd] = s;
				}
			}
		for (int i = 0; i < fns.size(); i++)
			fns[i]();
	}

	void print_help() {
		std::cerr << "usage:" << endl;
		std::cerr << "./search -c cache_file [-s init_cache] [-d 2,3...]" << endl;
		std::cerr << "parameter:" << endl;
		for (int i = 0; i < helps.size(); i++)
			std::cerr << helps[i] << endl;
		exit(1);
	}
} config;

struct Avg {
	LL total = 0, cnt = 0;
	void add(LL x) { total += x; cnt++; }
	double avg() { return (double)total / cnt; }
};

struct Logger {
	LL time;
	LL total_check_count = 0;
	LL total_fail_add_count = 0;
	// Avg shake;
	LL steps = 0;

	ofstream fout;
	void init() {
		time = gettime();
		string filename = config.log_filename;
		if (filename.empty())
			filename = "log_" + to_string(time) + ".txt";
		fout.open(filename, std::fstream::out);
		std::cerr << "writing log to " << filename << std::endl;
	}
	template<typename T>
	Logger& operator << (const T &t) {
		fout << t;
		std::cerr << t;
		return *this;
	}
	typedef Logger& (*LoggerManipulator)(Logger&);
	Logger& operator<<(LoggerManipulator manip) {
		return manip(*this);
	}
	static Logger& endl(Logger& logger) {
		std::cerr << std::endl;
		logger.fout << std::endl;
		return logger;
	}
	void log_pts(const vector<Point> &pt) {
		fout << pt.size() << ":";
		for(int i=1;i<=pt.size();i++) 
			fout<<"("<<pt[i-1].x<<","<<pt[i-1].y<<")";
	}
	void tic() {
		(*this) << "time:" << (gettime() - time);
		(*this) << " steps:" << steps;
		// (*this) << " avg_shake:" << shake.avg();
		(*this) << " check_count: " << total_check_count;
		(*this) << " fail_add_count: " << total_fail_add_count;
		(*this) << " @ " << log_beams() << Logger::endl;
	}
} logger;

bool check(const vector<Point> &pt, const Point &p){
	#pragma omp atomic
	logger.total_check_count++;
	return ahdoc::check(pt,p);
}

class Serializer {
public:
	std::fstream file;
public:
	Serializer() {}
	Serializer(const std::string &filename) : file(filename, std::ios::out | std::ofstream::binary) {}

	Serializer& operator <<(LL t) {
		file.write((char *)&t, sizeof(t));
		return *this;
	}
	template<typename T>
	Serializer& operator <<(std::vector<T> &vec) {
		(*this) << (LL)vec.size();
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
	Deserializer& operator >>(LL &t) {
		file.read((char *)&t, sizeof(t));
		return *this;
	}
	template<typename T>
	Deserializer& operator >>(std::vector<T> &vec) {
		LL n;
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

vector<Point> remove(vector<Point> pt, int D = 3, LL T = 100) {
	vector<Point> tmp;
	while (1) {
		tmp = sample(pt, pt.size() - D);
		if (ahdoc::check_no6hole(tmp)) {
			return tmp;
		}
	}
	return vector<Point>();
}

/* shake */

vector<int> PolarAngleStructure(vector<Point> pt_, int t) {
	thread_local vector<int> empty;
	thread_local vector<int> ret;
	thread_local vector<Point> ptsk;
	int n = pt_.size();
	ptsk.clear();
	ret.clear();
	for(int i=0;i<n;i++)
		if(i!=t){
			if(pt_[i]==pt_[t])
				return empty;
			ptsk.emplace_back(pt_[i].x - pt_[t].x , pt_[i].y - pt_[t].y);
		}
	--n;
	for(int i=0;i<n;i++)
		ret.push_back(i);
	sort(ret.begin(),ret.end(),[&](int i,int j){return ahdoc::atan2_less(ptsk[i],ptsk[j]);});
	
	// if (n >= 2 && !ahdoc::check_linear(ptsk.data(), n))
	// 	return empty;
	for(int i=1;i<n;i++)
		if(ptsk[ret[i-1]].x*ptsk[ret[i]].y==ptsk[ret[i]].x*ptsk[ret[i-1]].y)
			return empty;
	for(int i=0,j=1;i<n;i++){
		while(ptsk[ret[i]].x*ptsk[ret[j]].y-ptsk[ret[j]].x*ptsk[ret[i]].y>0) if (++j == n) j = 0;
		if(i!=j && ptsk[ret[i]].x*ptsk[ret[j]].y==ptsk[ret[j]].x*ptsk[ret[i]].y)
			return empty;
	}
	int st = 0;
	while (ret[st] != 0) st++;
	rotate(ret.begin(),ret.begin()+st,ret.end());

	for(int i=0,j=1;i<n;i++){
		while(ptsk[ret[i]].x*ptsk[ret[j]].y-ptsk[ret[j]].x*ptsk[ret[i]].y>0) if (++j == n) j = 0;
		ret.push_back(ret[j]);
	}
	
	return ret;
}

vector<Point> local_shake(vector<Point> pt, LL T = 0) {
	if(T<0) return pt;
	static const LL dx[4] = {-1, 0, 0, 1};
	static const LL dy[4] = {0, -1, 1, 0};
	vector<Point> tmp(pt.size());
	for(;T>=0;T--){
		tmp = pt;
		int i=randint(0,pt.size()-1);
		int chosen = randint(0,3);
		tmp[i].x += dx[chosen];
		tmp[i].y += dy[chosen];
		if(PolarAngleStructure(pt,i)!=PolarAngleStructure(tmp,i)) continue;
		//if (ahdoc::check_no6hole(tmp)) 
			pt=tmp;
		//else
		//	cerr<<"ERROR!\n";
	}
	return pt;
}

vector<Point> shake(vector<Point> pt, LL T = 0) {
	if(T<0) return pt;
	static const LL dx[4] = {-1, 0, 0, 1};
	static const LL dy[4] = {0, -1, 1, 0};
	vector<Point> tmp(pt.size());
	for(;T>=0;T--){
		tmp = pt;
		for(int i=0;i<pt.size()-1;i++){
			int chosen = randint(-4, 3);
			if(chosen>=0){
				tmp[i].x += (1<<T)*dx[chosen];
				tmp[i].y += (1<<T)*dy[chosen];
			}
		}
		if (ahdoc::check_no6hole(tmp))
			pt=tmp;
	}
	return pt;
}

vector<Point> expand(vector<Point> pt, LL R = 5000, LL T = 1000) {
	Point p;
	while (T--) {
		p.x = randint(-R, R);
		p.y = randint(-R, R);
		if (check(pt, p)) {
			pt.push_back(p);
			return pt;
		}
		//if(pt.size()>=24)
			// pt=shake(pt,1);
	}
	return vector<Point>();
}
vector<Point> expands_xx(vector<Point> pt, LL I = 1) {
	// repeat expand at most I times
	vector<Point> tmp;
	LL init = pt.size();
	while (I--) {
		LL R = config.r * pow(1.5874,(double)pt.size());
		//if(pt.size()>=25) R*=2;
		//if(pt.size()>=26) R*=3;
		//if(pt.size()>=27) R*=4;
		//if(pt.size()>=28) R*=5;
		//if(pt.size()>=29) R*=6;
		LL T = 1000;//(pt.size()<=20) ? (3LL*(LL)pt.size()*(LL)pt.size()) : 1000;
		tmp = expand(pt, R, T);
		if (tmp.empty())
			return vector<Point>();
		pt = tmp;
	}
	return pt;
}

struct BeamSearch {
	struct State {
		vector<Point> pt;
		int history = 0;
		Deserializer& Deserialize(Deserializer &ar) {
			return ar >> pt;
		}
		Serializer& Serialize(Serializer &ar) {
			return ar << pt;
		}
	};
	Serializer& Serialize(Serializer &ar) {
		return ar << width << largest << largest_pt << Q;
	}
	Deserializer& Deserialize(Deserializer &ar) {
		return ar >> width >> largest >> largest_pt >> Q;
	}
	vector<State> Q;
	LL head, width, largest, steps;
	vector<Point> largest_pt;
	std::string id;

	BeamSearch(LL width = 1000, const vector<Point> &pt = vector<Point>({{0, 0}}), std::string id = "") : width(width), id(id) {
		largest = 0;
		head = 0;
		Add(State({pt}));
	}

	vector<Point> pt_;
	State tmp;
	bool dfs(int i, int A) {
		if (A == 0) {
			vector<Point> pt = expands_xx(tmp.pt, pt_.size() - tmp.pt.size());
			if (!pt.empty()) {
				vector<Point> xxx = expands_xx(pt, 1);
				if (!xxx.empty()) {
					if (Add(State({xxx}))) return true;
				}
				else 
					Add(State({pt}));
			}
		}else{
			for (int j = i; j < pt_.size(); j++) {
				if (!check(tmp.pt, pt_[j]))
					continue;
				tmp.pt.push_back(pt_[j]);
				if (dfs(j + 1, A - 1))
					return true;
				tmp.pt.pop_back();
			}
		}
		return false;
	}
	
	int flatten_iter;
	vector<bool> flatten_pointers;
	void flatten_initialise_pointers(int n,int c){
		flatten_pointers.resize(n);
		for(int i=0;i<n;i++) flatten_pointers[i]=(i<c?true:false);
	}
	bool flatten_next_pointers(int n,int c){
		for(int i=n-1,cnt=0;i>=0;i--){
			if(flatten_pointers[i]) ++cnt;
			else{
				if(cnt==c) return false;
				for(int j=i-1;j>=0;j--)
					if(flatten_pointers[j]){
						flatten_pointers[j]=false;
						for(int k=j+1;k<n;k++)
							flatten_pointers[k]=(k-j<=cnt+1?true:false);
						return true;
					}
			}
		}
		return false;
	}
	
	bool Add(const State &stat) {
		// return true if stat.pt.size() > largest
		if (!stat.pt.empty()) {			
			if (largest < stat.pt.size()) {
				Q.clear();
				Q.push_back(stat);
				sort(Q[0].pt.begin(),Q[0].pt.end());
				head = 0;
				flatten_iter=0; flatten_initialise_pointers(Q[0].pt.size(),config.dels[flatten_iter]);
				
				if (stat.pt.size() >= 25) {
					#pragma omp critical
					{
						logger.log_pts(stat.pt);
						logger.fout << " by " << id << " @ " << steps;
						logger.fout << std::endl;
						// if (!largest_pt.empty()) {
						// 	logger.fout << "from: ";
						// 	logger.log_pts(largest_pt);
						// 	logger.fout << std::endl;
						// }
					}
				}
				if (stat.pt.size() >= 30)
				#pragma omp critical
				{
					logger << "Found point sets with " << stat.pt.size() << " points!!! checking..." << Logger::endl;
					if (!ahdoc::check_no6hole_bruteforce(stat.pt, true)) {
						std::cerr << "Error!" << endl;
						exit(1);
					}
					else {
						std::cerr << "Success!" << endl;
						exit(0);
					}
				}
				largest = stat.pt.size();
				largest_pt = stat.pt;
				return true;
			}
			else if (largest == stat.pt.size()) {
				if (Q.size() < width) {
					Q.push_back(stat);
				}
				else {
					Q[randint(0, Q.size() - 1)] = stat;
					#pragma omp atomic
					logger.total_fail_add_count++;
				}
			}
		}
		return false;
	}

	int pop() {
		int ret = head++;
		if (head == Q.size())
			head = 0;
		return ret;
	}

	void step() {
		steps++;
		if (largest < 16) {
			const vector<Point> &tmp = expands_xx(vector<Point>({ {0,0} }), 16);
			Add(State({tmp}));
		} else if (flatten_iter<config.dels.size()) {
			pt_ = Q[0].pt;
			tmp.pt.clear();
			for(int i=0;i<pt_.size();i++){
				if(!flatten_pointers[i])
					tmp.pt.push_back(pt_[i]);
			}
			if(ahdoc::check_no6hole(tmp.pt)){
				vector<Point> pt = expands_xx(tmp.pt,config.dels[flatten_iter]);
				if (!pt.empty()) {
					vector<Point> xxx = expands_xx(pt, 1);
					if (!xxx.empty()){
						Add(State({xxx}));
						return; // Exit the step() 
					}else
						Add(State({pt}));
				}
			}
			if(!flatten_next_pointers(Q[0].pt.size(),config.dels[flatten_iter])){
				++flatten_iter;
				if(flatten_iter<config.dels.size()) flatten_initialise_pointers(Q[0].pt.size(),config.dels[flatten_iter]);
			}
		}
		else {
			int sid = pop();
			Q[sid].pt = shake(Q[sid].pt, config.k);
			vector<Point> xxx = expands_xx(Q[sid].pt, 1);
			if (!xxx.empty())
				Add(State({xxx}));
			else {
				LL del_iters = config.dels.size();
				for (int iter = 0; iter < del_iters; iter++) {
					int c = config.dels[iter];
					vector<Point> tmp = expands_xx(remove(Q[sid].pt, c), c);
					if (!tmp.empty()) {
						vector<Point> xxx = expands_xx(tmp);
						if (!xxx.empty())
							Add(State({ xxx }));
						else
							Q[sid].pt = tmp;
						break;
					}
				}
			}
		}
	}
	void run() {
		while (1)
			step();
	}
};

vector<BeamSearch> ps;
string log_beams() {
	vector<LL> cnt;
	for (int i = 0; i < ps.size(); i++) {
		if (ps[i].largest >= cnt.size())
			cnt.resize(ps[i].largest + 1, 0);
		cnt[ps[i].largest]++;
	}
	string ret = "";
	for (int i = 0; i < cnt.size(); i++)
		if (cnt[i] > 0) {
			if (ret != "")
				ret += " ";
			ret += to_string(i) + ":" + to_string(cnt[i]);
		}
	return ret;
}

#define FOR(i, l, r) for(int i(l); i <(int)(r); i++)
int main(int argc, char *argv[]) {
	// std::ios::sync_with_stdio(false);
	config.parse(argc, argv);
	logger.init();
	logger << "Args: ";
	for (int i = 0; i < argc; i++)
		logger << string(argv[i]) << " ";
	logger << Logger::endl;
	logger << "Config: " << config.to_string() << Logger::endl;
	omp_set_num_threads(config.t);

	if (config.start_cache_filename != "") {
		logger << "start loading from " << config.start_cache_filename << Logger::endl;
		Deserializer(config.start_cache_filename) >> ps;
		logger << "progress loaded from " << config.start_cache_filename << " @ " << log_beams() << Logger::endl;
		if (config.b > 1) {
			LL m = 0;
			for (int i = 0; i < ps.size(); i++)
				if (ps[i].largest >= config.b)
					ps[m++] = ps[i];
			ps.resize(m);
			logger << "progress trimed to @ " << log_beams() << Logger::endl;
			
			vector<LL> cnt;
			for (int i = 0; i < ps.size(); i++) {
				if (ps[i].head >= cnt.size())
					cnt.resize(ps[i].head + 1, 0);
				cnt[ps[i].head]++;
			}
			string ret = "";
			for (int i = 0; i < cnt.size(); i++)
				if (cnt[i] > 0) {
					if (ret != "")
						ret += " ";
					ret += to_string(i) + ":" + to_string(cnt[i]);
				}
			logger << "progress heads @ " << ret << Logger::endl;
		}
	}
	else if (config.start_log_filename != "") {
		logger << "start loading from " << config.start_log_filename << Logger::endl;
		ifstream fin(config.start_log_filename);
		string line;
		std::vector<std::vector<Point>> pts;
		// auto peek = [&](const std::string &s, int i, int len) {
		// 	if (i + len <= s.length())
		// 		return s.substr(i, len);
		// 	else
		// 		return std::string();
		// };
		auto read_pts = [&](const std::string &s) {
			istringstream ssin(line);
			int n; 
			LL x, y;
			char c;
			ssin >> n;
			vector<Point> pt;
			for (int i = 0; i < n; i++) {
				while (c != '(') ssin >> c;
				ssin >> x >> c >> y >> c;
				pt.push_back({x, y});
			}
			return pt;
		};
		while (getline(fin, line)) {
			if (line.length() < 3) continue;
			if (line[2] == ':' && isdigit(line[0]) && isdigit(line[1])) {
				std::vector<Point> pt = read_pts(line);
				if (config.lt != -1 && pt.size() != config.lt) continue;
				if (!ahdoc::check_no6hole(pt)) continue;
				pts.push_back(pt);
			}
		}
		std::reverse(pts.begin(), pts.end());
		for (int i = 0; i < pts.size(); i++) {
			ps.push_back(BeamSearch(config.w, pts[i], "ps_" + std::to_string(i)));
			if (ps.size() >= config.n)
				break;
		}
		logger << "progress loaded from " << config.start_log_filename << " @ " << log_beams() << Logger::endl;
	}
	else {
		ps.resize(config.n, BeamSearch(config.w));
	}
	for (int i = 0; i < ps.size(); i++)
		ps[i].id = "ps_" + std::to_string(i);
	
	LL last_cache_time = gettime(), last_log_time = gettime(), cache_id = 0;
	LL interval = 10 * 1000; // quick log in the first a fwe tics
	while (1) {
		logger.steps++;
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < ps.size(); i++) {
			ps[i].step();
			if (omp_get_thread_num() == 0 && (gettime() - last_log_time) >= interval) 
			#pragma omp critical
			{
				logger.tic();
				last_log_time = gettime();
				interval = min(interval + 10 * 1000, config.log_interval * 1000);
			}
		}
		if (config.cache_filename != "") {
			if ((gettime() - last_cache_time) >= config.cache_interval * 1000) {
				string filename = config.cache_filename + "_" + to_string(cache_id % 2 + 1) + ".binary";
				logger << "start saving to " << filename << " @ " << log_beams() << Logger::endl;
				Serializer(filename) << ps;
				logger << "progress saved to " << filename << Logger::endl;
				cache_id++;
				last_cache_time = gettime();
			}
		}
	}
}

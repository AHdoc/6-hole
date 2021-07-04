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
		int s1 = (((p1.y > 0) << 1) - 1) << (p1.x < 0);
		int s2 = (((p2.y > 0) << 1) - 1) << (p2.x < 0);
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
	
	bool find6hole(const vector<Point> &_pt, const Point &p){
		n=_pt.size();
		for(int i=0;i<n;i++){
			pt[i].x=_pt[i].x-p.x;
			pt[i].y=_pt[i].y-p.y;
		}
		sort(pt,pt+n,atan2_less);
		if(solve()) return true;
		for(int i=0;i<n;i++){
			pt[i].x=p.x-_pt[i].x;
			pt[i].y=p.y-_pt[i].y;
		}
		sort(pt,pt+n,atan2_less);
		if(solve()) return true;
		return false;
	}
	
	bool check(const vector<Point> &pt, const Point &p){
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

	bool check_no6hole(vector<Point> pt){ 
		// Check whether or not the given set of points contains no 6-hole
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
	LL n = 100, w = 1000, t = 7, b = 0, k = 1;
	LL cache_interval = 3600, log_interval = 300;
	string log_filename = "", cache_filename = "", start_filename = "", start_log_filename = "";
	vector<int> dels = {3,3};

	string to_string() {
		string ret;
		ret += "n=" + std::to_string(n) + ", ";
		ret += "w=" + std::to_string(w) + ", ";
		ret += "t=" + std::to_string(t) + ", ";
		ret += "k=" + std::to_string(k) + ", ";
		ret += "log_filename=" + log_filename + ", ";
		ret += "cache_filename=" + cache_filename + ", ";
		ret += "start_filename=" + start_filename + ", ";
		ret += "start_log_filename=" + start_log_filename + ", ";
		ret += "b=" + std::to_string(b) + ", ";
		ret += "dels=";
		for (int i = 0; i < dels.size(); i++)
			ret += (i?",":"") + std::to_string(dels[i]);
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
				if (cmd == "o")
					log_filename = peek(i + 1);
				else if (cmd == "c")
					cache_filename = peek(i + 1);
				else if (cmd == "l")
					start_log_filename = peek(i + 1);
				else if (cmd == "s")
					start_filename = peek(i + 1);
				else if (cmd == "b")
					b = stoll(peek(i + 1));
				else if (cmd == "k")
					k = stoll(peek(i + 1));
				else if (cmd == "n")
					n = stoll(peek(i + 1));
				else if (cmd == "t")
					t = stoll(peek(i + 1));
				else if (cmd == "w")
					w = stoll(peek(i + 1));
				else if (cmd == "d")
					dels = parse_vec(peek(i + 1));
				else if (cmd == "i")
					cache_interval = stoll(peek(i + 1));
				else if (cmd == "t")
					log_interval = stoll(peek(i + 1));
				else if (cmd == "h")
					print_help();
			}
	}

	void print_help() {
		std::cerr << "usage:" << endl;
		std::cerr << "./search -c cache_file [-s init_cache] [-d 2,3...]" << endl;
		std::cerr << "parameter:" << endl;
		std::cerr << "    -o string: output file [default: log_#.txt]" << endl;
		std::cerr << "    -c string: cache file, no cache if empty [default: \"\"]" << endl;
		std::cerr << "    -l string: init log file, no init loading if empty [default: \"\"]" << endl;
		std::cerr << "    -s string: init cache file, no init loading if empty [default: \"\"]" << endl;
		std::cerr << "    -n integer: number of beams [default: 3000]" << endl;
		std::cerr << "    -b integer: number of least loaded pts [default: 1]" << endl;
		std::cerr << "    -k integer: number of shakes [default: 10]" << endl;
		std::cerr << "    -t integer: number of threads [default: 7]" << endl;
		std::cerr << "    -w integer: size of beams [default: 1000]" << endl;
		std::cerr << "    -d list: number points to delete [default: 3,3]" << endl;
		std::cerr << "    -i integer: number of seconds between cache [default: 3600]" << endl;
		std::cerr << "    -l integer: max number of seconds between logs [default: 300]" << endl;
		std::cerr << "    -h: print help message and quit" << endl;
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
		fout<<std::endl;
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

vector<Point> shake(vector<Point> pt, LL R = 1, LL T = 1000) {
	// vector<Point> pt_ = pt;
	static const LL dx[4] = {-1, 0, 0, 1};
	static const LL dy[4] = {0, -1, 1, 0};
	vector<Point> tmp(pt.size());
	for (int iter = 0; iter < T; iter++) {
		int i = randint(0, pt.size()-1);
		tmp = pt;
		int chosen = randint(0, 3);
		tmp[i].x += dx[chosen];
		tmp[i].y += dy[chosen];
		if (ahdoc::check_no6hole(tmp))
			pt=tmp;
	}
	// LL dist = 0;
	// for (int i = 0; i < pt.size(); i++) {
	// 	dist += abs(pt[i].x - pt_[i].x) + abs(pt[i].y - pt_[i].y);
	// }
	// #pragma omp critical 
	// logger.shake.add(dist);
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
	}
	return vector<Point>();
}
vector<Point> expands_xx(vector<Point> pt, LL I = 1) {
	// repeat expand at most I times
	vector<Point> tmp;
	LL init = pt.size();
	while (I--) {
		LL R = 100 * pt.size();
		LL T = (pt.size()<=20) ? (3LL*(LL)pt.size()*(LL)pt.size()) : (1000 + max((LL)pt.size() - 20, 0LL)*max((LL)pt.size() - 20, 0LL) * 1000);
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
	LL head, width, largest;
	vector<Point> largest_pt;

	BeamSearch(LL width = 1000, const vector<Point> &pt = vector<Point>({{0, 0}})) : width(width) {
		largest = 0;
		head = 0;
		Add(State({pt}));
	}

	void Add(const State &stat) {
		if (!stat.pt.empty()) {
			if (largest < stat.pt.size()) {
				largest = stat.pt.size();
				largest_pt = stat.pt;
				Q.clear();
				Q.push_back(stat);
				head = 0;
				if (stat.pt.size() >= 25) {
					#pragma omp critical
					logger.log_pts(stat.pt);
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
	}

	int pop() {
		int ret = head++;
		if (head == Q.size())
			head = 0;
		return ret;
	}

	void step() {
		if (largest < 21) {
			const vector<Point> &tmp = expands_xx(vector<Point>({ {0,0} }), 21);
			Add(State({tmp}));
		} else {
			int sid = pop();
			vector<Point> pt = Q[sid].pt;
			// add
			LL add_iters = 1;
			for (int iter = 0; iter < add_iters; iter++) {
				const vector<Point> &tmp = expands_xx(pt, 1);
				Add(State({tmp}));
			}
			// del
			LL del_iters = 0;
			if (pt.size() >= 21) {
				del_iters = 2;
			}
			del_iters = min(del_iters, (LL)config.dels.size());
			for (int iter = 0; iter < del_iters; iter++) {
				int c = config.dels[iter];
				const vector<Point> &tmp = expands_xx(shake(remove(pt, c),1,config.k), c);
				const vector<Point> &xxx = expands_xx(tmp, 1);
				if (xxx.empty())
					Add(State({tmp}));
				else
					Add(State({xxx}));
			}
			//shake itself
			if (pt.size() >= 22 && largest == pt.size()) {
				Q[sid] = State({ shake(pt, 1, config.k) });
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

	if (config.start_filename != "") {
		logger << "start loading from " << config.start_filename << Logger::endl;
		Deserializer(config.start_filename) >> ps;
		logger << "progress loaded from " << config.start_filename << " @ " << log_beams() << Logger::endl;
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
		while (getline(fin, line)) {
			if (line.length() < 3) continue;
			if (line[2] == ':' && isdigit(line[0]) && isdigit(line[1])) {
				istringstream ssin(line);
				int n; 
				LL x, y;
				char c;
				ssin >> n >> c;
				if(n<26) continue;
				vector<Point> pt;
				for (int i = 0; i < n; i++) {
					ssin >> c >> x >> c >> y >> c;
					pt.push_back({x, y});
				}
				ps.push_back(BeamSearch(config.w, pt));
			}
		}
		logger << "progress loaded from " << config.start_log_filename << " @ " << log_beams() << Logger::endl;
	}
	else {
		ps.resize(config.n, BeamSearch(config.w));
	}
	
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

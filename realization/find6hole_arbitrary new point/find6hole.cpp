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

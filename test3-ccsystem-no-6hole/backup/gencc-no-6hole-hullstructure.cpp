/*
DIMACS format:
	http://www.satcompetition.org/2009/format-benchmarks2009.html
CC system:
	https://en.wikipedia.org/wiki/CC_system

At least 30 points are needed:
	there exists a set of 29 points in general position with no empty convex hexagon.
*/

#include<iostream>
#include<fstream>
#include<cstdio>
#include<cstring>
#include<vector>
#include<set>
#include<algorithm>

#include<assert.h> 

using namespace std;

const int MAXN=50;

int nbhullstructure;
int nbvar,nbclauses,nbliterals;

int lvl[MAXN+1],known[MAXN+1][MAXN+1][MAXN+1],known2[MAXN+1][MAXN+1][MAXN+1][MAXN+1];
int idx[MAXN+1][MAXN+1][MAXN+1];
int idx2[MAXN+1][MAXN+1][MAXN+1][MAXN+1];
set<vector<int>> clauses;

pair<int,int> oppo(pair<int,int> c){
	return make_pair(-c.first,-c.second);
}

void new_clauses(vector<pair<int,int>> c){
	vector<int> d;
	for(auto pii:c){
		if(pii.second==1) return;
		if(pii.second==-1) continue;
		d.push_back(pii.first);
	}
	
	assert(!d.empty());
	
	sort(d.begin(),d.end());
	if(clauses.find(d)==clauses.end()){
		++nbclauses;
		nbliterals+=d.size();
		clauses.insert(d);
	}
}

void define_var_and(pair<int,int> nv,vector<pair<int,int>> cond){
	vector<pair<int,int>> tmp;
	for(auto x:cond)
		tmp.push_back(oppo(x));
	tmp.push_back(nv);
	new_clauses(tmp);
	
	for(auto x:cond)
		new_clauses({oppo(nv),x});
} 

pair<int,int> get(int i,int j,int k){
	assert(i!=j && i!=k && j!=k);
	if(i<j && i<k){
		if(j<k) return make_pair(idx[i][j][k],known[i][j][k]);
		else return make_pair(-idx[i][k][j],-known[i][k][j]);
	}else
		return get(j,k,i);
}

void find_6hole(int nbpt,int st,int ed,int k,vector<int> pts){
	if(k==0){
		vector<pair<int,int>> tmp;
		int n=pts.size();
		for(int i=0;i<n;i++){
			int previ=(i+n-1)%n;
			int j=(i+1)%n;
			for(int u=0;u<n;u++)
				if(u!=i && u!=j && u!=previ)
					tmp.push_back(oppo(get(pts[i],pts[j],pts[u])));
		}
		for(int i=1;i<=nbpt;i++){
			if(find(pts.begin(),pts.end(),i)==pts.end()){
				tmp.push_back(make_pair(idx2[pts[0]][pts[2]][pts[3]][i],known2[pts[0]][pts[2]][pts[3]][i]));
				tmp.push_back(make_pair(idx2[pts[0]][pts[3]][pts[4]][i],known2[pts[0]][pts[3]][pts[4]][i]));
			}
		}
		new_clauses(tmp);
		return;
	}
	vector<int> pts2;
	for(int i=st;i<=ed;i++)
		if(find(pts.begin(),pts.end(),i)==pts.end()){
			pts2=pts;
			pts2.push_back(i);
			find_6hole(nbpt,st,ed,k-1,pts2);
		}
}

void cc_system(int n){
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				if(known[i][j][k]==0){
					++nbvar;
					idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=nbvar;
					idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-nbvar;
				}
			}
	// Interiority: If tqr and ptr and pqt, then pqr.
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int t=1;t<=n;t++){
					set<int> pts={p,q,r,t};
					if(pts.size()!=4) continue;
					new_clauses({oppo(get(t,q,r)),oppo(get(p,t,r)),oppo(get(p,q,t)),get(p,q,r)});
				}
	// Transitivity: If tsp and tsq and tsr, and tpq and tqr, then tpr.
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++)
					for(int t=1;t<=n;t++){
						set<int> pts={p,q,r,s,t};
						if(pts.size()!=5) continue;
						new_clauses({oppo(get(t,s,p)),oppo(get(t,s,q)),oppo(get(t,s,r)),oppo(get(t,p,q)),oppo(get(t,q,r)),get(t,p,r)});
					}
}

void var_pt_inside_triangle(int n){
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++)
					idx2[p][q][r][s]=-1;
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					if(known2[p][q][r][s]!=0) continue;
					if(idx2[p][q][r][s]!=-1) continue;
					idx2[p][q][r][s]=idx2[q][r][p][s]=idx2[r][p][q][s]=++nbvar;
					define_var_and(make_pair(idx2[p][q][r][s],known2[p][q][r][s]),{get(p,q,r),get(p,q,s),get(q,r,s),get(r,p,s)});
				}
}

void mk_no6hole_givenhullstructure(int n,vector<int> hull){
	nbvar=0;
	nbclauses=0;
	nbliterals=0;
	clauses.clear();
	
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			for(int k=1;k<=n;k++)
				known[i][j][k]=0;
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++)
					known2[p][q][r][s]=0;
	for(int i=1,i2,j=0;j<hull.size();i=i2,j++){
		i2=i+hull[j];
		vector<int> layerpts;
		for(int k=i;k<i2;k++){
			lvl[k]=j;
			layerpts.push_back(k);
		}
		int tot=layerpts.size();
		if(tot>=3){
			for(int u=0;u<tot;u++)
				for(int v=(u+1)%tot;(v+1)%tot!=u;v=(v+1)%tot)
					for(int w=(v+1)%tot;w!=u;w=(w+1)%tot){
						int uu=layerpts[u];
						int vv=layerpts[v];
						int ww=layerpts[w];
						known[uu][vv][ww]=1;
						known[uu][ww][vv]=-1;
					}
			for(int u=0;u<tot;u++){
				int v=(u+1)%tot;
				int p=layerpts[u],q=layerpts[v];
				for(int r=i2;r<=n;r++){
					known[p][q][r]=known[q][r][p]=known[r][p][q]=1;
					known[p][r][q]=known[q][p][r]=known[r][q][p]=-1;
				}
			}
		}
	}
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					if(lvl[s]<=lvl[p] && lvl[s]<=lvl[q] && lvl[s]<=lvl[r])
						known2[p][q][r][s]=-1;
				}
	
	cc_system(n);
	var_pt_inside_triangle(n);
	
	// Restriction: No 6-hole.
	for(int p=1;p<=n;p++)
		find_6hole(n,p+1,n,5,{p});
	
	/* 
	string s=to_string(n)+"pts-no-6hole-";
	for(int x:hull)
		s.push_back(char('0'+x));
	s+=".sat";
	ofstream fout(s);
	
	fout<<"c\n";
	fout<<"p cnf "<<nbvar<<" "<<clauses.size()<<"\n"; 
	for(auto elm:clauses){
		for(int x:elm)
			fout<<x<<" ";
		fout<<"0\n";
	}
	fout.close();
	*/
	
	cerr<<n<<" size ("<<nbvar<<","<<nbclauses<<","<<nbliterals<<")\n";
}

void mk_no6hole_enumerate_hullstructure(int n,int sumx,int bigx,vector<int> hull){
	if(sumx==n){
		if(hull[hull.size()-1]>=6) return; // [1]
		
		++ nbhullstructure;
		for(int x:hull)
			cout<<x<<" ";
		cout<<"\n";
		
		mk_no6hole_givenhullstructure(n,hull);
	}else if(sumx+1==n){
		if(hull.size()>=1 && hull[hull.size()-1]<=7){ // [2]
			vector<int> newhull=hull;
			newhull.push_back(1);
			mk_no6hole_enumerate_hullstructure(n,n,max(1,bigx),newhull);
		}
	}else if(sumx+2==n){
		if(hull.size()>=1 && hull[hull.size()-1]<=6){ // [3]
			vector<int> newhull=hull;
			newhull.push_back(2);
			mk_no6hole_enumerate_hullstructure(n,n,max(2,bigx),newhull);
		}
	}else{
		for(int i=3;i<9 && i+sumx<=n;i++){
			if(i==3 && hull.size()>=1 && hull[hull.size()-1]>7) continue; // [4]
			
			vector<int> newhull=hull;
			newhull.push_back(i);
			mk_no6hole_enumerate_hullstructure(n,sumx+i,max(i,bigx),newhull);
		}
	}
}

int main(){
	int n;
	cin>>n;
	nbhullstructure=0;
	mk_no6hole_enumerate_hullstructure(n,0,0,{});
	cerr<<"nbhullstructure = "<<nbhullstructure<<"\n";
} 


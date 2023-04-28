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
#include<unordered_set>
#include<algorithm>

#include<assert.h> 

using namespace std;

struct VectorHash {
    size_t operator()(const vector<int>& v) const {
        hash<int> hasher;
        size_t seed = 0;
        for (int i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

const int MAXN=30;

int nbhullstructures;
int nbvar,nbclauses,nbliterals;

int idx[MAXN+1][MAXN+1][MAXN+1];
int idx2[MAXN+1][MAXN+1][MAXN+1][MAXN+1]; 
int idx3[MAXN+1][MAXN+1][MAXN+1];
int idx4[MAXN+1][MAXN+1][MAXN+1][MAXN+1]; 
int lvl[MAXN+1];

unordered_set<vector<int>,VectorHash> clauses;

int new_var(){
	++nbvar;
	return nbvar;
}

void new_clauses(vector<int> d){
	sort(d.begin(),d.end());
	if(clauses.find(d)==clauses.end()){
		++nbclauses;
		nbliterals+=d.size();
		clauses.insert(d);
	}
}

void define_var_and(int nv,vector<int> cond){
	vector<int> tmp;
	for(int x:cond) tmp.push_back(-x);
	tmp.push_back(nv);
	new_clauses(tmp);
	
	for(int x:cond)	new_clauses({-nv,x});
}

int get(int i,int j,int k){
	assert(i!=j && i!=k && j!=k);
	if(i<j && i<k){
		if(j<k) return idx[i][j][k];
		else return -idx[i][k][j];
	}else
		return get(j,k,i);
}

void new_known(int x){
	new_clauses({x});
}

void new_known(int i,int j,int k){
	new_known(idx[i][j][k]);
}

void cc_system(int n,vector<int> hull){
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				int x=new_var();
				idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=x;
				idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-x;
			}
/* known information */
	for(int i=1,i2,j=0;j<hull.size();i=i2,j++){ 
		if(j>=1 && hull[j-1]==3 && hull[j]>=2){ // Symmetry Breaking for Triangles
			int p=i-3,q=i-2,x=i,y=i+1;
			new_known(p,y,x);
			new_known(q,x,y);
		} 
		if(j>=1 && hull[j-1]==4){ // Symmetry Breaking for Quadrilaterals
			int p=i-4,q=i-3,r=i-2,s=i-1;
			new_known(q,s,i);
			new_known(r,p,i);
		}
		if(j>=1 && hull[j-1]==5){ // Symmetry Breaking for Pentagons
			int p=i-5,q=i-4,r=i-3,s=i-2,t=i-1;
			int x1=new_var(),x6=new_var(),x11=new_var();
			new_known(q,s,i);
			new_known(s,p,i);
			new_known(r,t,i);
			new_clauses({get(q,t,i),get(p,r,i)}); 
		} 
		if(j>=1 && hull[j-1]==6){ // Symmetry Breaking for Hexagons
			int p=i-6,q=i-5,r=i-4,s=i-3,t=i-2,u=i-1;
			int x1=new_var(),x6=new_var(),x11=new_var();
			new_known(q,t,i);
			new_known(s,p,i);
		}
		
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
						new_known(uu,vv,ww);
					}
			for(int u=0;u<tot;u++){
				int v=(u+1)%tot;
				int p=layerpts[u],q=layerpts[v];
				for(int r=i2;r<=n;r++)
					new_known(p,q,r);
			}
		}
	}
/* ----------------- */
	// Interiority: If tqr and ptr and pqt, then pqr.
	for(int p=1;p<=n;p++)
		for(int q=p+1;q<=n;q++)
			for(int r=p+1;r<=n;r++)
				for(int t=1;t<=n;t++){
					set<int> pts={p,q,r,t};
					if(pts.size()!=4) continue;
					new_clauses({-get(t,q,r),-get(p,t,r),-get(p,q,t),get(p,q,r)});
				}
	// Transitivity: If tsp and tsq and tsr, and tpq and tqr, then tpr.
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++)
					for(int t=1;t<=n;t++){
						set<int> pts={p,q,r,s,t};
						if(pts.size()!=5) continue;
						new_clauses({-get(t,s,p),-get(t,s,q),-get(t,s,r),-get(t,p,q),-get(t,q,r),get(t,p,r)});
					}
}

void var_pt_inside_triangle(int n,vector<int> hull){
	for(int p=1;p<=n;p++)
		for(int q=p+1;q<=n;q++)
			for(int r=p+1;r<=n;r++){
				vector<int> empty_cond;
				for(int s=1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					
					int x=new_var();
					if(lvl[s]<=lvl[p] && lvl[s]<=lvl[q] && lvl[s]<=lvl[r]) new_known(-x);
					
					idx2[p][q][r][s]=idx2[q][r][p][s]=idx2[r][p][q][s]=x;
					define_var_and(x,{get(p,q,s),get(q,r,s),get(r,p,s)});
					
					empty_cond.push_back(-x);
				}
				int x_empty=new_var();
				idx3[p][q][r]=idx3[q][r][p]=idx3[r][p][q]=x_empty;
				define_var_and(x_empty,empty_cond); 
			}
}

void var_4pt_region_empty(int n,vector<int> hull){
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					
					vector<int> empty_cond;
					for(int t=1;t<=n;t++){
						if(pts.find(t)!=pts.end()) continue;
						
						int x=new_var(); 
						define_var_and(x,{get(p,q,t),get(r,s,t),get(p,s,t)});
						empty_cond.push_back(-x);
					}
					int y=new_var();
					idx4[p][q][r][s]=y;
					define_var_and(y,empty_cond);
				}
}

void mk_no6hole_givenhullstructure(int n,vector<int> hull){
	nbvar=0;
	nbclauses=0;
	nbliterals=0;
	clauses.clear();
	
	cc_system(n,hull);
	var_pt_inside_triangle(n,hull);
	var_4pt_region_empty(n,hull);
	
	// Restriction: No 6-hull.
	for(int p=1;p<=n;p++)
		for(int q=p+1;q<=n;q++)
			for(int r=p+1;r<=n;r++)
				for(int s=p+1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					
					int empty_r1=idx3[p][q][r],empty_r2=idx3[p][r][s];
					int left_empty=idx4[p][q][r][s],right_empty=idx4[r][s][p][q];
					new_clauses({-get(p,q,r),-get(q,r,s),-get(r,s,p),-get(s,p,q),-empty_r1,-empty_r2,left_empty,right_empty});
				}
	
	string s=to_string(n)+"pt-no-6hole-";
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
	
	cerr<<n<<" ("<<nbvar<<","<<nbclauses<<","<<nbliterals<<")\n";
}

/* 8=0*
   81=81*
   82=82*
   83*
   84*
   85=850*
   852=852*
   853*
*/ 

vector<vector<int>> banned={
	//6<=n<=15 
	{6},
	{7},
	{8},
	{7,2},{8,1},
	{8,2},
	{6,3,2},{7,3,1},{8,3},
	{7,3,2},{8,3,1},{8,4},
	{7,3,3},{8,3,2},{8,4,1},{8,5},
	{6,3,3,2},{6,3,4,1},{6,3,5},{7,3,3,1},{7,3,4},{7,4,3},{8,3,3},{8,4,2},
	{6,4,3,2},{7,3,3,2},{7,3,4,1},{7,3,5},{7,4,3,1},{8,3,3,1},{8,3,4},{8,4,3},{8,5,2},
	//n=16
	{5,3,3,3,2},{6,3,3,3,1},{6,3,3,4},{7,3,3,3},{7,3,4,2},{7,3,5,1},{7,4,3,2},{7,4,5},{7,5,3,1},
	{8,3,3,2},{8,3,4,1},{8,3,5},{8,4,3,1},{8,4,4},{8,5,3},{8,7,1},
	//n=17 and hull[0]=8
	{8,3,3,3},{8,3,4,2},{8,3,5,1},{8,4,3,2},{8,4,4,1},{8,4,5},{8,5,3,1},{8,5,4}
};

bool banned_hull(vector<int> hull){
	int n=hull.size();
	for(auto pat:banned){
		int m=pat.size();
		if(n>=m){
			bool suf=true;
			for(int i=0;i<m;i++)
				if(hull[n-m+i]!=pat[i])
					suf=false;
			if(suf)
				return true;
		}
	}
	return false;
}

void mk_no6hole_enumerate_hullstructure(int n,int sumx,vector<int> hull){
	if(sumx==n){
		++nbhullstructures;
		for(int x:hull)
			cout<<x<<" ";
		cout<<"\n";
		
		if(!banned_hull(hull) && hull[0]==8)
			mk_no6hole_givenhullstructure(n,hull);
	}else if(sumx+1==n){
		vector<int> newhull=hull;
		newhull.push_back(1);
		mk_no6hole_enumerate_hullstructure(n,n,newhull);
	}else if(sumx+2==n){
		vector<int> newhull=hull;
		newhull.push_back(2);
		mk_no6hole_enumerate_hullstructure(n,n,newhull);
	}else{
		for(int i=3;i<9 && i+sumx<=n;i++){
			vector<int> newhull=hull;
			newhull.push_back(i);
			mk_no6hole_enumerate_hullstructure(n,sumx+i,newhull);
		}
	}
}

int main(){
	mk_no6hole_givenhullstructure(8+5+3,{8,5,3});
	mk_no6hole_givenhullstructure(8+5+3+1,{8,5,3,1});
	mk_no6hole_givenhullstructure(8+5+3+2,{8,5,3,2});
	mk_no6hole_givenhullstructure(8+5+3+3,{8,5,3,3});
	mk_no6hole_givenhullstructure(8+5+3+4,{8,5,3,4});
	mk_no6hole_givenhullstructure(8+5+3+5,{8,5,3,5});
	
	//for(int n=17;n<=17;n++){
	//	nbhullstructures=0;
	//	mk_no6hole_enumerate_hullstructure(n,0,{});
	//}
} 


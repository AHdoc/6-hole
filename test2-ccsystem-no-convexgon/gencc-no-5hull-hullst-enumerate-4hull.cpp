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

const int MAXN=50;

int nbhullstructures;
int nbvar,nbclauses,nbliterals;

int idx[MAXN+1][MAXN+1][MAXN+1],known[MAXN+1][MAXN+1][MAXN+1];
int idx2[MAXN+1][MAXN+1][MAXN+1][MAXN+1]; 
int lvl[MAXN+1];

unordered_set<vector<int>,VectorHash> clauses;

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

void cc_system(int n,vector<int> hull){
/* known information */
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			for(int k=1;k<=n;k++)
				known[i][j][k]=0;
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
/* ----------------- */
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				++nbvar;
				idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=nbvar;
				idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-nbvar;
				
				if(known[i][j][k]==1)
					new_clauses({nbvar});
				if(known[i][j][k]==-1)
					new_clauses({-nbvar});
			}
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
						
						define_var_and(++nbvar,{get(p,q,t),get(r,s,t),get(p,s,t)});
						empty_cond.push_back(-nbvar);
					}
					idx2[p][q][r][s]=++nbvar;
					define_var_and(idx2[p][q][r][s],empty_cond);
				}
}

void mk_no5hull_givenhullstructure(int n,vector<int> hull){
	nbvar=0;
	nbclauses=0;
	nbliterals=0;
	clauses.clear();
	
	cc_system(n,hull);
	var_4pt_region_empty(n,hull);
	
	// Restriction: No 5-hull.
	for(int p=1;p<=n;p++)
		for(int q=p+1;q<=n;q++)
			for(int r=p+1;r<=n;r++)
				for(int s=p+1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					
					int left_empty=idx2[p][q][r][s];
					new_clauses({-get(p,q,r),-get(q,r,s),-get(r,s,p),-get(s,p,q),left_empty});
				}
	
	string s=to_string(n)+"pts-no-5hull-";
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

void mk_no5hull_enumerate_hullstructure(int n,int sumx,vector<int> hull){
	if(sumx==n){
		++nbhullstructures;
		for(int x:hull)
			cout<<x<<" ";
		cout<<"\n";
		
		mk_no5hull_givenhullstructure(n,hull);
	}else if(sumx+1==n){
		vector<int> newhull=hull;
		newhull.push_back(1);
		mk_no5hull_enumerate_hullstructure(n,n,newhull);
	}else if(sumx+2==n){
		vector<int> newhull=hull;
		newhull.push_back(2);
		mk_no5hull_enumerate_hullstructure(n,n,newhull);
	}else{
		for(int i=3;i<5 && i+sumx<=n;i++){
			vector<int> newhull=hull;
			newhull.push_back(i);
			mk_no5hull_enumerate_hullstructure(n,sumx+i,newhull);
		}
	}
}

int main(){
	for(int n=8;n<=9;n++){
		nbhullstructures=0;
		mk_no5hull_enumerate_hullstructure(n,0,{});
	}
} 


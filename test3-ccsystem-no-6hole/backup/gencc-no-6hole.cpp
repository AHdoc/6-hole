/*
DIMACS format:
	http://www.satcompetition.org/2009/format-benchmarks2009.html
CC system:
	https://en.wikipedia.org/wiki/CC_system

At least 30 points are needed:
	there exists a set of 29 points in general position with no empty convex hexagon.
*/

#include<iostream>
#include<cstdio>
#include<cstring>
#include<vector>
#include<set>
#include<algorithm>

#include<assert.h> 

using namespace std;

const int MAXN=50;

int nbvar,nbclauses,nbliterals;

int idx[MAXN+1][MAXN+1][MAXN+1];
int idx2[MAXN+1][MAXN+1][MAXN+1][MAXN+1];
vector<vector<int>> clauses;

void new_clauses(vector<int> d){
	++nbclauses;
	nbliterals+=d.size();
	//clauses.push_back(d);
}

int get(int i,int j,int k){
	assert(i!=j && i!=k && j!=k);
	if(i<j && i<k){
		if(j<k) return idx[i][j][k];
		else return -idx[i][k][j];
	}else
		return get(j,k,i);
}

void find_6hole(int nbpt,int st,int ed,int k,vector<int> pts){
	if(k==0){
		vector<int> tmp;
		int n=pts.size();
		for(int i=0;i<n;i++){
			int previ=(i+n-1)%n;
			int j=(i+1)%n;
			for(int u=0;u<n;u++)
				if(u!=i && u!=j && u!=previ)
					tmp.push_back(-get(pts[i],pts[j],pts[u]));
		}
		for(int i=1;i<=nbpt;i++){
			if(find(pts.begin(),pts.end(),i)==pts.end()){
				tmp.push_back(idx2[pts[0]][pts[1]][pts[2]][i]);
				tmp.push_back(idx2[pts[0]][pts[2]][pts[3]][i]);
				tmp.push_back(idx2[pts[0]][pts[3]][pts[4]][i]);
				tmp.push_back(idx2[pts[0]][pts[4]][pts[5]][i]);
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
				++nbvar;
				idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=nbvar;
				idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-nbvar;
			}
	// Interiority: If tqr and ptr and pqt, then pqr.
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
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
					if(idx2[p][q][r][s]!=-1) continue;
					idx2[p][q][r][s]=idx2[q][r][p][s]=idx2[r][p][q][s]=++nbvar;
					new_clauses({-get(p,q,r),-get(p,q,s),-get(q,r,s),-get(r,p,s),idx2[p][q][r][s]});
					new_clauses({-idx2[p][q][r][s],get(p,q,r)});
					new_clauses({-idx2[p][q][r][s],get(p,q,s)});
					new_clauses({-idx2[p][q][r][s],get(q,r,s)});
					new_clauses({-idx2[p][q][r][s],get(r,p,s)});
				}
}

void mk_no6hole(int n){
	nbvar=0;
	nbclauses=0;
	nbliterals=0;
	clauses.clear();
	
	cc_system(n);
	var_pt_inside_triangle(n);
	
	// Restriction: No 6-hole.
	for(int p=1;p<=n;p++)
		find_6hole(n,p+1,n,5,{p});
	
	//cout<<"c\n";
	//cout<<"p cnf "<<nbvar<<" "<<clauses.size()<<"\n"; 
	//for(auto elm:clauses){
	//	for(int x:elm)
	//		cout<<x<<" ";
	//	cout<<"0\n";
	//}
	
	cerr<<"nbclauses  = "<<nbclauses<<"\n";
	cerr<<"nbliterals = "<<nbliterals<<"\n";
}

int main(){
	int n;
	cin>>n;
	mk_no6hole(n);
} 


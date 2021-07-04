/*
DIMACS format:
	http://www.satcompetition.org/2009/format-benchmarks2009.html
CC system:
	https://en.wikipedia.org/wiki/CC_system
*/

#include<iostream>
#include<cstdio>
#include<cstring>
#include<vector>
#include<set>
#include<algorithm>
using namespace std;

const int MAXN=100;

int idx[MAXN+1][MAXN+1][MAXN+1];
vector<vector<int>> clauses;

int get(int i,int j,int k){
	if(i<j && i<k){
		if(j<k) return idx[i][j][k];
		else return -idx[i][k][j];
	}else
		return get(j,k,i);
}

void find_gon(int st,int ed,int k,vector<int> pts){
	if(k==0){
		vector<int> tmp;
		int n=pts.size();
		if(n==5){
			tmp.push_back(-get(pts[0],pts[1],pts[2]));
			tmp.push_back(-get(pts[0],pts[1],pts[3]));
			tmp.push_back(-get(pts[0],pts[1],pts[4]));
			tmp.push_back(-get(pts[1],pts[2],pts[3]));
			tmp.push_back(-get(pts[1],pts[2],pts[4]));
			tmp.push_back(-get(pts[2],pts[3],pts[4]));
			tmp.push_back(-get(pts[2],pts[3],pts[0]));
			tmp.push_back(-get(pts[3],pts[4],pts[0]));
		}else{
			for(int i=0;i<n;i++){
				int j=(i+1)%n;
				for(int u=0;u<n;u++)
					if(u!=i && u!=j)
						tmp.push_back(-get(pts[i],pts[j],pts[u]));
			}
		}
		clauses.push_back(tmp);
		return;
	}
	vector<int> pts2;
	for(int i=st;i<=ed;i++)
		if(find(pts.begin(),pts.end(),i)==pts.end()){
			pts2=pts;
			pts2.push_back(i);
			find_gon(st,ed,k-1,pts2);
		}
}

void mk(int n,int len){
	int nbvar=0;
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				++nbvar;
				idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=nbvar;
				idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-nbvar;
			}
	
	clauses.clear();
	// Interiority: If tqr and ptr and pqt, then pqr.
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int t=1;t<=n;t++){
					set<int> pts={p,q,r,t};
					if(pts.size()!=4) continue;
					clauses.push_back({-get(t,q,r),-get(p,t,r),-get(p,q,t),get(p,q,r)});
				}
	// Transitivity: If tsp and tsq and tsr, and tpq and tqr, then tpr.
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++)
					for(int t=1;t<=n;t++){
						set<int> pts={p,q,r,s,t};
						if(pts.size()!=5) continue;
						clauses.push_back({-get(t,s,p),-get(t,s,q),-get(t,s,r),-get(t,p,q),-get(t,q,r),get(t,p,r)});
					}
	
	// Restriction: No len-gon.
	for(int p=1;p<=n;p++)
		find_gon(p+1,n,len-1,{p});
	
	cout<<"c\n";
	cout<<"p cnf "<<nbvar<<" "<<clauses.size()<<"\n";
	for(auto elm:clauses){
		for(int x:elm)
			cout<<x<<" ";
		cout<<"0\n";
	}
}

int main(){
	int n,m;
	cin>>n>>m;
	mk(n,m);
} 

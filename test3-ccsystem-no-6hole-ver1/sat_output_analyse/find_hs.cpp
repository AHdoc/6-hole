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

const int MAXN=30;

int nbvar;

/* Variables */
int idx[MAXN+1][MAXN+1][MAXN+1]; // idx[p][q][r]: The triangle pqr is oriented
int state[MAXN*MAXN*MAXN+1];

int new_var(){ // Set a new variable
	++nbvar;
	return nbvar;
}

void cc_system(int n){
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				int x=new_var();
				idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=x;
				idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-x;
			}
}

void find_hs(int n){
	nbvar=0;
	cc_system(n);
	
	cerr<<"nbvar = "<<nbvar<<"\n";
	cerr<<"idx[3][2][4] = "<<idx[3][2][4]<<"\n";
	
	for(int i=1;i<=nbvar;i++){
		cin>>state[i];
		assert(abs(state[i])==i);
	}

	set<int> pt;
	for(int i=1;i<=n;i++) pt.insert(i);
	
	while(pt.size()>=3){
		set<int> hull;
		for(int i:pt){
			for(int j:pt){
				if(i==j) continue;
				bool ck=true;
				for(int k:pt){
					if(i==k || j==k) continue;
					assert(idx[i][j][k]!=0);
					if((idx[i][j][k]>0 && state[idx[i][j][k]]<0)||(idx[i][j][k]<0 && state[-idx[i][j][k]]>0)){
						ck=false;
						break;
					}
				}
				if(ck){
					cout<<i<<"->"<<j<<"\n";
					hull.insert(i);
					hull.insert(j);
				}
			}
		}
		cout<<"hull.size() = "<<hull.size()<<"\n";
		if(hull.size()==0) exit(1); 
		for(int x:hull){
			pt.erase(pt.find(x));
		}
	}
	if(!pt.empty()) cout<<"hull.size() = "<<pt.size()<<"\n";
}

int main(int argc, char *argv[]){
	freopen("output_sat_4466710_extra_1.txt","r",stdin);
	find_hs(29);
}


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

int nbvar,nbclauses,nbliterals;

/* Variables */
int idx[MAXN+1][MAXN+1][MAXN+1]; // idx[p][q][r]: The triangle pqr is oriented
int idx2[MAXN+1][MAXN+1][MAXN+1];
int idx3[MAXN+1][MAXN+1][MAXN+1][MAXN+1]; 

int lvl[MAXN+1];

unordered_set<vector<int>,VectorHash> clauses;

int new_var(){ // Set a new variable
	++nbvar;
	return nbvar;
}

void new_clauses(vector<int> d){ // Give a new clause
	sort(d.begin(),d.end());
	if(clauses.find(d)==clauses.end()){
		++nbclauses;
		nbliterals+=d.size();
		clauses.insert(d);
	}
}

void define_var_and(int nv,vector<int> cond){ // Define the new variable nv as (x1 & x2 & ... & xk) with cond={x1,...,xk}
	vector<int> tmp;
	for(int x:cond) tmp.push_back(-x);
	tmp.push_back(nv);
	new_clauses(tmp);
	
	for(int x:cond)	new_clauses({-nv,x});
}

void new_known(int x){ // Claim that the variable x is known where x may be negative
	new_clauses({x});
}

void new_known(int i,int j,int k){
	new_known(idx[i][j][k]);
}

void cc_system(int n){
	for(int i=1;i<=n;i++)
	for(int j=i+1;j<=n;j++)
	for(int k=j+1;k<=n;k++){
		int x=new_var();
		idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=x;
		idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-x;
	}
	// Interiority: If tqr and ptr and pqt, then pqr.
	for(int p=1;p<=n;p++)
	for(int q=1;q<=n;q++)
	for(int r=1;r<=n;r++)
	for(int t=1;t<=n;t++){
		set<int> pts={p,q,r,t};
		if(pts.size()!=4) continue;
		new_clauses({-idx[t][q][r],-idx[p][t][r],-idx[p][q][t],idx[p][q][r]});
	}
	// Transitivity: If tsp and tsq and tsr, and tpq and tqr, then tpr.
	for(int p=1;p<=n;p++)
	for(int q=1;q<=n;q++)
	for(int r=1;r<=n;r++)
	for(int s=1;s<=n;s++)
	for(int t=1;t<=n;t++){
		set<int> pts={p,q,r,s,t};
		if(pts.size()!=5) continue;
		new_clauses({-idx[t][s][p],-idx[t][s][q],-idx[t][s][r],-idx[t][p][q],-idx[t][q][r],idx[t][p][r]});
	}
}

// idx2[p][q][r]: the triangle pqr has no point inside
void var_pt_inside_triangle(int n){
	for(int p=1;p<=n;p++)
	for(int q=1;q<=n;q++)
	for(int r=1;r<=n;r++){
		set<int> pts1={p,q,r};
		if(pts1.size()!=3) continue;

		vector<int> empty_cond;
		for(int s=1;s<=n;s++){
			set<int> pts2={p,q,r,s};
			if(pts2.size()!=4) continue;
			
			int x=new_var(); // x: s is inside triangle pqr
			define_var_and(x,{idx[p][q][s],idx[q][r][s],idx[r][p][s]});
			
			empty_cond.push_back(-x);
		}

		int x_empty=new_var();
		idx2[p][q][r]=x_empty;
		define_var_and(x_empty,empty_cond); 
	}
}

// idx3[p][q][r][s]: A special 4-point region is empty
/*   -------s------r
     xxxxxxx|
	 xxxxxxx|
	 xxxxxxx|
	 -------p------q 
*/ 
void var_4pt_region_empty(int n){
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
			define_var_and(x,{idx[p][q][t],idx[r][s][t],idx[p][s][t]});
			empty_cond.push_back(-x);
		}

		int y=new_var();
		idx3[p][q][r][s]=y;
		define_var_and(y,empty_cond);
	}
}

/*
1:   2:   3:   4:   5:   6:   7:   8:
 +A+  A-+  +-A  +A+  +A+  A-+  +-A  
 B C  | C  B |  | C  B |  | |  | |   A
 +D+  +D+  +D+  B-+  +-C  +-C  B-+
*/

const int num_points[9]={-1,4,3,3,3,3,2,2,1};

void clauses_from_hull(int n,vector<int> hull){
	int outside_n=2,iA,iB,iC,iD; // 1 : plus x inf; 2 : plus y inf.
	for(int x:hull){
		if(x==1){
			iA=outside_n+1, iB=outside_n+2, iC=outside_n+3, iD=outside_n+4;
		}else if(x==2){
			iA=outside_n+1, iB=outside_n+1, iC=outside_n+2, iD=outside_n+3;
		}else if(x==3){
			iA=outside_n+1, iB=outside_n+2, iC=outside_n+1, iD=outside_n+3;
		}else if(x==4){
			iA=outside_n+1, iB=outside_n+2, iC=outside_n+3, iD=outside_n+2;			
		}else if(x==5){
			iA=outside_n+1, iB=outside_n+2, iC=outside_n+3, iD=outside_n+3;			
		}else if(x==6){
			iA=outside_n+1, iB=outside_n+1, iC=outside_n+2, iD=outside_n+2;
		}else if(x==7){
			iA=outside_n+1, iB=outside_n+2, iC=outside_n+1, iD=outside_n+2;
		}else if(x==8){
			iA=outside_n+1, iB=outside_n+1, iC=outside_n+1, iD=outside_n+1;			
		}
		for(int i=outside_n+1;i<=n;i++) if(i!=iA) new_known(idx[1][iA][i]);
		for(int i=outside_n+1;i<=n;i++) if(i!=iB) new_known(idx[2][iB][i]);
		for(int i=outside_n+1;i<=n;i++) if(i!=iC) new_known(idx[iC][2][i]);
		for(int i=outside_n+1;i<=n;i++) if(i!=iD) new_known(idx[iD][1][i]);
		outside_n+=num_points[x];
	}
}

void mk_no6hole_given_recthull(vector<int> hull){
	int n=2; for(int x:hull) n+=num_points[x]; // Total number of points
	nbvar=0; nbclauses=0; nbliterals=0; clauses.clear(); // Initialise the SAT problem
	
	cc_system(n);
	clauses_from_hull(n,hull);
	var_pt_inside_triangle(n);
	var_4pt_region_empty(n);
	
	// Restriction: No 6-hole.
	for(int p=1;p<=n;p++)
		for(int q=p+1;q<=n;q++)
			for(int r=p+1;r<=n;r++)
				for(int s=p+1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					
					int empty_r1=idx2[p][q][r],empty_r2=idx2[p][r][s];
					int left_empty=idx3[p][q][r][s],right_empty=idx3[r][s][p][q];
					new_clauses({-idx[p][q][r],-idx[q][r][s],-idx[r][s][p],-idx[s][p][q],-empty_r1,-empty_r2,left_empty,right_empty});
				}
	
	// Output SAT problem
	string s="no-6hole-";
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

int main(){
	mk_no6hole_given_recthull({1,1,1,8});
}

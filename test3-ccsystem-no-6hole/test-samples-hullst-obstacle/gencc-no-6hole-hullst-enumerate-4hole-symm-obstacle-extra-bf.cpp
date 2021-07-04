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
int idx3[MAXN+1][MAXN+1][MAXN+1];
int idx4[MAXN+1][MAXN+1][MAXN+1][MAXN+1]; 
int idx5[MAXN+1][MAXN+1][MAXN+1][MAXN+1];
int idx6[MAXN+1][MAXN+1][MAXN+1][MAXN+1]; 

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

void define_var_or(int nv,vector<int> cond){ // Define the new variable nv as (x1 | x2 | ... | xk) with cond={x1,...,xk}
	vector<int> tmp;
	for(int x:cond) tmp.push_back(x);
	tmp.push_back(-nv);
	new_clauses(tmp);
	
	for(int x:cond)	new_clauses({nv,-x});
}

void new_known(int x){ // Claim that the variable x is known where x may be negative
	new_clauses({x});
}

void new_known(int i,int j,int k){
	new_known(idx[i][j][k]);
}

void cc_system(int n,vector<int> hull,int kn){
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				int x=new_var();
				idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=x;
				idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-x;
			}
/* known information */
	for(int i=1;i<=n;i++)
		lvl[i]=-1;
	for(int i=kn+1,i2,j=0;j<hull.size();i=i2,j++){
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
	// 1 ~ kn are outside
	//for(int p=1;p<=kn;p++){
	//	vector<int> tmp;
	//	for(int i=kn+1;i<=kn+hull[0];i++){
	//		int j=(i==kn+hull[0]?kn+1:i+1);
	//		tmp.push_back(idx[j][i][p]);
	//	}
	//	new_clauses(tmp);
	//}
	
/* ----------------- */
	// Interiority: If tqr and ptr and pqt, then pqr.
	for(int p=1;p<=n;p++)
		for(int q=p+1;q<=n;q++)
			for(int r=p+1;r<=n;r++)
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
// idx5[p][q][r][s]: pq and rs (along the counterclockwise order prqs) intersect
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					int x=new_var();
					idx5[p][q][r][s]=x;
					define_var_and(x,{idx[p][q][s],idx[q][p][r],idx[r][s][p],idx[s][r][q]});
				}
// idx6[p][q][r][s]: pq and the left extension of rs intersect
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					int x=new_var();
					idx6[p][q][r][s]=x;
					define_var_and(x,{idx[p][r][q],idx[p][s][q],idx[r][s][q],idx[s][r][p]});
				}
}

int var_triangle_not_intersect_obstable(int n,int p,int q,int r,vector<int> hull){
	int x=new_var();
	if(lvl[p]==hull.size()-1 && lvl[q]==hull.size()-1 && lvl[r]==hull.size()-1){
		new_known(-x);
		return x;
	}
	vector<int> cond;
	int lt=n-hull[hull.size()-1]+1,rt=n;
	for(int i1=lt;i1<=rt;i1++){
		int i2=(i1<rt?i1+1:lt);
		if((set<int>{p,q,i1,i2}).size()==4) cond.push_back(-idx5[p][q][i1][i2]);
		if((set<int>{q,r,i1,i2}).size()==4) cond.push_back(-idx5[q][r][i1][i2]);
		if((set<int>{r,p,i1,i2}).size()==4) cond.push_back(-idx5[r][p][i1][i2]);
		if((set<int>{p,q,i1,i2}).size()==4) cond.push_back(-idx5[p][q][i2][i1]);
		if((set<int>{q,r,i1,i2}).size()==4) cond.push_back(-idx5[q][r][i2][i1]);
		if((set<int>{r,p,i1,i2}).size()==4) cond.push_back(-idx5[r][p][i2][i1]);
	}
	define_var_and(x,cond);
	return x;
}
// idx3[p][q][r]   : the triangle pqr has no point inside nor intersection with the obstacle, i,e. is totally empty
void var_pt_inside_triangle(int n,vector<int> hull){
	for(int p=1;p<=n;p++)
		for(int q=p+1;q<=n;q++)
			for(int r=p+1;r<=n;r++){
				if(q==r) continue; // To guarantee that p, q and r are distinct

				vector<int> empty_cond;
				for(int s=1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue; // To guarantee that p, q, r and s are distinct
					
					int x=new_var(); // x: the point [s] is inside the triangle pqr
					define_var_and(x,{idx[p][q][s],idx[q][r][s],idx[r][p][s]});
					if(lvl[s]!=-1 && lvl[s]<=lvl[p] && lvl[s]<=lvl[q] && lvl[s]<=lvl[r]) new_known(-x); // Pruning
					
					empty_cond.push_back(-x);
				}
				int x_not_intersect_obstacle=var_triangle_not_intersect_obstable(n,p,q,r,hull);
				empty_cond.push_back(x_not_intersect_obstacle);
				
				int x_empty=new_var();
				idx3[p][q][r]=idx3[q][r][p]=idx3[r][p][q]=x_empty;
				define_var_and(x_empty,empty_cond); 
			}
}

// idx4[p][q][r][s]: A special 4-point region is NOT empty (in any case) 
/*   -------s------r
     xxxxxxx|
	 xxxxxxx|
	 xxxxxxx|
	 -------p------q 
*/ 
void var_4pt_region_nonempty(int n,vector<int> hull){
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;  // To guarantee that p, q, r and s are distinct
					
					vector<int> nonempty_cond;
					for(int t=1;t<=n;t++){
						if(pts.find(t)!=pts.end()) continue;  // To guarantee that p, q, r, s and t are distinct
						
						int x=new_var(); // x: the point [t] is inside the 4-point region pqrs 
						define_var_and(x,{idx[p][q][t],idx[r][s][t],idx[p][s][t]});
						
						nonempty_cond.push_back(x);
					}
					int y=new_var();
					idx4[p][q][r][s]=y;
					define_var_or(y,nonempty_cond);
				}
}

void mk_no6hole_givenhullstructure_obstacle(vector<int> hull,int kn){
	int n=kn;
	for(int x:hull) n+=x;
	
	nbvar=0;
	nbclauses=0;
	nbliterals=0;
	clauses.clear();
	
	cc_system(n,hull,kn);
	var_pt_inside_triangle(n,hull);
	var_4pt_region_nonempty(n,hull);
	
	// Restriction: No 6-hole.
	for(int p=1;p<=n;p++)
		for(int q=p+1;q<=n;q++)
			for(int r=p+1;r<=n;r++)
				for(int s=p+1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					
					int empty_r1=idx3[p][q][r],empty_r2=idx3[p][r][s];
					int left_nonempty=idx4[p][q][r][s],right_nonempty=idx4[r][s][p][q];
					new_clauses({-idx[p][q][r],-idx[q][r][s],-idx[r][s][p],-idx[s][p][q],-empty_r1,-empty_r2,-left_nonempty,-right_nonempty});
				}
	
	string s="no-6hole-";
	for(int x:hull)
		s.push_back(char('0'+x));
	s+="-obstacle";
	if(kn>0)
		s+="-"+to_string(kn)+"points-extra"; 
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

int main(int argc, char *argv[]){
	if(argv[1][0]=='-'){
		if(argv[1][1]=='h'){
			printf("No one can help you.\n");
		}
	}else{
		if(argc>=2){
			vector<int> hullst;
			
			int len=strlen(argv[1]),n=0;
			for(int i=0;i<len;i++){
				int x=argv[1][i]-'0';
				hullst.push_back(x);
				n+=x;
			}
			
			int k=0;
			if(argc>=3)
				k=atoi(argv[2]);
			n+=k;
			
			printf("n = %d.\n", n);
			if(k>0) printf("k = %d points extra.\n", k);
			printf("Given hull structure = {%d",hullst[0]);
			for(int i=1;i<hullst.size();i++) printf(",%d",hullst[i]);
			printf("}.\n");
			
			mk_no6hole_givenhullstructure_obstacle(hullst,k);
		}
	}
}


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

ifstream fin;
ofstream fout;

/* Variables */
int idx[MAXN+1][MAXN+1][MAXN+1]; // idx[p][q][r]: The triangle pqr is oriented
unordered_set<vector<int>,VectorHash> clauses;

void new_clauses(vector<int> d){ // Give a new clause
	sort(d.begin(),d.end());
	clauses.insert(d);
}

void new_known(int x){ // Claim that the variable x is known where x may be negative
	new_clauses({x});
}

void new_known(int i,int j,int k){
	new_known(idx[i][j][k]);
}

void cc_system(int k,int m,vector<int> secondhull,string c,int s){
	int n=k+m;
	int cnt=0;
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				++cnt;
				idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=cnt;
				idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-cnt;
			}
	
/* known information */
	// restrict the smallest hull in the second hull structure; we assume that n is the right most one in the smallest hull
	int sz=secondhull[secondhull.size()-1]; // points on the smallest hull are labelled n-sz+1 to n
	if(c=="covered"){
		for(int i=n-sz+1;i<=n;i++){
			int j=(i<n?i+1:n-sz+1);
			new_known(8,i,j);
		}
		int s2=(s<7?s+1:1);
		new_known(8,s,n);
		new_known(s2,8,n);
		for(int i=1;i<s;i++){
			int j=i+1;
			for(int k=n-sz+1;k<n;k++)
				new_clauses({idx[i][8][k],idx[8][j][k]});
		}
		for(int k=n-sz+1;k<n;k++)
			new_clauses({idx[s][8][k],idx[8][n][k]});
	}else if(c=="isolated"){
		int i=s; 
		int j=(i<7?i+1:1);
		new_known(8,i,n);
		new_known(j,8,n);
		for(int k=n-sz+1;k<n;k++)
			new_known(8,n,k);
	}
	
	// make points outside apart based on the second hull structure
	for(int i=k+1,i2,j=0;j<secondhull.size();i=i2,j++){ 
		if(j>=1 && secondhull[j-1]==3 && secondhull[j]>=2){ // Symmetry Breaking for Triangles
			int p=i-3,q=i-2,x=i,y=i+1;
			new_known(p,y,x);
			new_known(q,x,y);
		} 
		if(j>=1 && secondhull[j-1]==4){ // Symmetry Breaking for Quadrilaterals
			int p=i-4,q=i-3,r=i-2,s=i-1;
			new_known(q,s,i);
			new_known(r,p,i);
		}
		if(j>=1 && secondhull[j-1]==5){ // Symmetry Breaking for Pentagons
			int p=i-5,q=i-4,r=i-3,s=i-2,t=i-1;
			new_known(q,s,i);
			new_known(s,p,i);
			new_known(r,t,i);
			new_clauses({idx[q][t][i],idx[p][r][i]}); 
		} 
		if(j>=1 && secondhull[j-1]==6){ // Symmetry Breaking for Hexagons
			int p=i-6,q=i-5,r=i-4,s=i-3,t=i-2,u=i-1;
			new_known(q,t,i);
			new_known(s,p,i);
		}
		if(j>=1 && secondhull[j-1]==7){ // Symmetry Breaking for Heptagon
			int p=i-7,q=i-6,r=i-5,s=i-4,t=i-3,u=i-2,v=i-1;
			new_known(q,t,i);
			new_known(s,v,i);
			new_known(u,p,i);
			new_clauses({idx[q][v][i],idx[p][r][i]});
			new_clauses({idx[q][u][i],idx[p][s][i]});
		}
		
		i2=i+secondhull[j];
		vector<int> layerpts2;
		for(int k=i;k<i2;k++) layerpts2.push_back(k);
		int tot=layerpts2.size();
		if(tot>=3){
			for(int u=0;u<tot;u++)
				for(int v=(u+1)%tot;(v+1)%tot!=u;v=(v+1)%tot)
					for(int w=(v+1)%tot;w!=u;w=(w+1)%tot){
						int uu=layerpts2[u],vv=layerpts2[v],ww=layerpts2[w];
						new_known(idx[uu][vv][ww]);
					}
			for(int u=0;u<tot;u++){
				int v=(u+1)%tot;
				int p=layerpts2[u],q=layerpts2[v];
				for(int r=i2;r<=n;r++)
					new_known(idx[p][q][r]);
			}
		}
	}
}

void additional_secondhullst(int k,int m,string inputfile,vector<int> secondhullst,string c,int s){
	int n=k+m;
	
	string outputfile=inputfile+"-secondhs{";
	for(int i=0;i<secondhullst.size();i++) outputfile.push_back('0'+secondhullst[i]);
	outputfile+="-"+c+"-"+to_string(s)+"}.sat";
	inputfile+=".sat";
	
	cerr<<" input file = "<< inputfile<<"\n";
	cerr<<"output file = "<<outputfile<<"\n";
	
	fin.open(inputfile);
	fout.open(outputfile);
	
	string tmp1,tmp2;
	int nbvar,nbclause;
	fin>>tmp1;
	if(tmp1=="c") fin>>tmp1;
	fin>>tmp2>>nbvar>>nbclause;
	
	clauses.clear();
	cc_system(k,m,secondhullst,c,s);
	
	fout<<"p cnf "<<nbvar<<" "<<nbclause+clauses.size()<<"\n"; 
	for(auto elm:clauses){
		for(int x:elm)
			fout<<x<<" ";
		fout<<"0\n";
	}
	
	
	int tmp3;
	while(fin>>tmp3){
		while(tmp3!=0){
			fout<<tmp3<<" ";
			fin>>tmp3;
		}
		fout<<"0\n";
	}
	
	fin.close();
	fout.close();
}

int main(int argc, char *argv[]){
	ios_base::sync_with_stdio(false); 
	
	if(argv[1][0]=='-'){
		if(argv[1][1]=='h'){
			fprintf(stderr,
"  -h             print this short list of common options\n"
"usage: k m filename hull-structure c s\n"
"                 there are $k points in the middle labelled 1 to k, and $m points outside, labelled k+1 to k+m;\n"
"                 the clauses of the known part in the middle are given by the SAT file of name \"$filename\"+\".sat\";\n"
"                 the second hull structure for those points outside is given by $hull-structure, which is a string with characters '1'~'8'.\n"
"                 c=``covered'' or c=``isolated'':\n"
"                    [covered] : point 8 is covered by the smallest hull in the given hull structure;\n" 
"                                then s is the first non-empty region of the 7-partition.\n"
"                    [isolated]: point 8 is away from the smallest hull in the given hull structure;\n" 
"                                then s is the region of the 7-partition containing the right most point in the smallest hull.\n"
"                 ATTENTION: the size of the given hull-structure must be equal to $m.\n"
			);
		}
	}else{
		if(argc==7){
			int k=atoi(argv[1]);
			int m=atoi(argv[2]);
			string filename=argv[3];
			vector<int> secondhullst;
			
			int len=strlen(argv[4]),cnt=0;
			for(int i=0;i<len;i++){
				int x=argv[4][i]-'0';
				secondhullst.push_back(x);
				cnt+=x;
			}
			
			string c=argv[5];
			int s=atoi(argv[6]);
			
			if(cnt!=m){
				printf("The total size of the second hull structure = %d, but m = %d.\n",cnt,m);
				exit(1);
			}
			
			if(c=="covered" || c=="isolated"){
				if(1<=s && s<=7);
				else{
					printf("Invalid s.\n");
					exit(1);
				}
			}else{
				printf("Invalid c.\n");
				exit(1);
			} 
			
			printf("k = %d, m = %d, n = %d.\n",k,m,k+m);
			printf("second hull structure = {%d",secondhullst[0]);
			for(int i=1;i<secondhullst.size();i++) printf(",%d",secondhullst[i]);
			printf("}.\n");
			
			additional_secondhullst(k,m,filename,secondhullst,c,s);
		}
	}
} 




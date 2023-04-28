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

const int MAXN=30;

ifstream fin;
ofstream fout;

/* Variables */
int idx[MAXN+1][MAXN+1][MAXN+1]; // idx[p][q][r]: The triangle pqr is oriented

void new_clauses(vector<int> d){ // Give a new clause
	sort(d.begin(),d.end());
	
	for(int x:d)
		fout<<x<<" ";
	fout<<"0\n";
}

void new_known(int x){ // Claim that the variable x is known where x may be negative
	new_clauses({x});
}

void cc_system(int k,int m,vector<pair<bool,int>> partition){
	int n=k+m;
	int cnt=0;
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				++cnt;
				idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=cnt;
				idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-cnt;
			}
	
	vector<int> layerpts;
	for(int k=1;k<=7;k++) layerpts.push_back(k);
	
	for(int i=0,j,u=k+1,v;i<partition.size();i++,u=v){
		v=u+abs(partition[i].second);
		j=(i+1)%7;
		if(partition[i].first){
			for(int w=k+1;w<=n;w++){
				if(u<=w && w<v){
					new_known(idx[8][layerpts[i]][w]);
					new_known(-idx[8][layerpts[j]][w]);
				}else{
					new_clauses({-idx[8][layerpts[i]][w],idx[8][layerpts[j]][w]});
				}
			}
		}
	}
}

void additional_7parts(int k,int m,string inputfile,vector<pair<bool,int>> partition){
	string outputfile=inputfile+"(";
	for(int i=0;i<partition.size();i++){
		if(i!=0) outputfile+=",";
		if(!partition[i].first) outputfile+="x";
		outputfile+=to_string(partition[i].second);
	} 
	outputfile+=").sat";
	inputfile+=".sat";
	
	cerr<<" input file = "<< inputfile<<"\n";
	cerr<<"output file = "<<outputfile<<"\n";
	
	fin.open(inputfile);
	fout.open(outputfile);
	
	string tmp1,tmp2;
	int nbvar,nbclause,actual_nbclause;
	
	fin>>tmp1>>tmp2>>nbvar>>nbclause;
	actual_nbclause=nbclause;
	for(int i=0;i<7;i++)
		if(partition[i].first)
			actual_nbclause+=m+partition[i].second;
	
	fout<<tmp1<<" "<<tmp2<<" "<<nbvar<<" "<<actual_nbclause<<"\n";
	
	cc_system(k,m,partition);
	
	for(int i=1;i<=nbclause;i++){
		int tmp3;
		while(fin>>tmp3 && tmp3!=0){
			fout<<tmp3<<" ";
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
"usage: k m filename a b c d e f g\n"
"                 there are $k points in the middle labelled 1 to k, and $m points outside, labelled k+1 to k+m;\n"
"                 the clauses of the known part in the middle are given by the SAT file of name \"$filename\"+\".sat\";\n"
"                 the 7-partition is given by {a,...,g}, where only those without leading 'x' are asked to be added into the SAT file.\n"
"                 You may omit the last some but not all of the seven integers.\n" 
			);
		}
	}else{
		if(argc>11 || argc<=4){
			fprintf(stderr,"Error!\n");
		}else{
			int k=atoi(argv[1]);
			int m=atoi(argv[2]);
			string filenames=argv[3];
			
			vector<pair<bool,int>> partition;
			int cnt=0;
			for(int i=4;i<argc;i++){
				int x=(argv[i][0]=='x'?atoi(argv[i]+1):atoi(argv[i]));
				cnt+=abs(x);
				partition.push_back(make_pair(argv[i][0]!='x',x));
			}
			if(cnt>m || (partition.size()==7 && cnt<m)){
				printf("Error!\n");
			}else{
				while(partition.size()!=7)
					partition.push_back(make_pair(false,0));
				if(cnt<m)
					partition[partition.size()-1].second+=m-cnt;
				additional_7parts(k,m,filenames,partition);
			}
		}
	}
} 




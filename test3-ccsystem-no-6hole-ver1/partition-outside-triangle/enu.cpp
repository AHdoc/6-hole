#include<iostream>
#include<cstdio>
#include<cstring>
#include<set>
#include<vector>
#include<algorithm>
using namespace std;

void solve(int o){
	set<vector<int>> pat;
	for(int a=0;a<=o;a++)
	for(int b=0;b<=o;b++)
	for(int c=0;c<=o;c++)
	for(int d=0;d<=o;d++)
	for(int e=0;e<=o;e++)
	for(int f=0;f<=o;f++){
		if(a+b+c+d+e+f==o){
			if(f+a+b==0 || f+a+b==o || f+a+b>3 || a>2) continue;
			if(b+c+d==0 || b+c+d==o || b+c+d>3 || c>2) continue;
			if(d+e+f==0 || d+e+f==o || d+e+f>3 || e>2) continue;
			
			if(f+a+b==1 || b+c+d==1 || d+e+f==1) continue;
			if(f+a+b==2 || b+c+d==2 || d+e+f==2) continue;
			
			vector<vector<int>> symm;
			symm.push_back({a,b,c,d,e,f});
			symm.push_back({c,d,e,f,a,b});
			symm.push_back({e,f,a,b,c,d});
			symm.push_back({a,f,e,d,c,b});
			symm.push_back({e,d,c,b,a,f});
			symm.push_back({c,b,a,f,e,d});
			sort(symm.begin(),symm.end());
			pat.insert(symm[0]);
		}
	}
	for(auto x:pat){
		for(auto y:x) cout<<y<<" ";
		cout<<"\n";
	}
}

int main(){
	solve(4);
}

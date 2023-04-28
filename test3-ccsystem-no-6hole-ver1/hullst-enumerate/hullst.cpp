#include<iostream>
#include<cstdio>
#include<cstring>
using namespace std;

int cnt;

void check(int n,string s,bool smallest_hulls,bool hulls_with_obstacle,int startswith_l,int startswith_r){
	if(startswith_l<=s[0]-'0' && s[0]-'0'<=startswith_r){
		if(hulls_with_obstacle){
			if(s[s.size()-1]=='1' || s[s.size()-1]=='2');
			else cout<<++cnt<<" "<<n<<" "<<s<<"\n";
		}
		if(smallest_hulls) cout<<++cnt<<" "<<n<<" "<<s<<"0\n";
	}
}

void dfs(int tot,int n,string s,bool smallest_hulls,bool hulls_with_obstacle,int startswith_l,int startswith_r){
	if(n<=2){
		if(n>0) s.push_back('0'+n);
		check(tot,s,smallest_hulls,hulls_with_obstacle,startswith_l,startswith_r);
	}else{
		//for(int i=3;i<=8 && i<=n;i++){
		for(int i=min(8,n);i>=3;i--){
			string t=s;
			t.push_back('0'+i);
			dfs(tot,n-i,t,smallest_hulls,hulls_with_obstacle,startswith_l,startswith_r);
		}
	}
}

int main(){
	/*
	freopen("n_26_startswith(4).txt","w",stdout);
	int st_l=4,st_r=4;
	cnt=0;
	for(int n=26;n<=26;n++){
		//if(n<=27) dfs(n,n,"",true,true,st_l,st_r);
		//else 
		dfs(n,n,"",true,false,st_l,st_r);
	}
	fclose(stdout);
	*/
	
	for(int n=1;n<=30;n++){
		freopen(("n_"+to_string(n)+"_smallest_hulls.txt").c_str(),"w",stdout);
		cnt=0;
		dfs(n,n,"",true,false,3,8);
		fclose(stdout);
	} 
}


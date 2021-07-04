#include<iostream>
#include<cstdio>
#include<cstring>
#include<vector>
#include<set>
#include<algorithm>
using namespace std;

int main(){
	int n;
	cin>>n;
	set<vector<int>> S;
	for(int a=0;a<=n;a++)
		for(int b=0;b<=n;b++)
			for(int c=0;c<=n;c++)
				for(int d=0;d<=n;d++)
					for(int e=0;e<=n;e++)
						for(int f=0;f<=n;f++)
							for(int g=0;g<=n;g++){
								if(a+b+c+d+e+f+g!=n) continue;
								vector<int> ff={a,b,c,d,e,f,g};
								vector<vector<int>> hh;
								for(int i=0;i<7;i++){
									vector<int> ff2;
									for(int j=0;j<7;j++) ff2.push_back(ff[(i+j)%7]);
									hh.push_back(ff2);
								}
								for(int i=0;i<7;i++){
									vector<int> ff2;
									for(int j=0;j<7;j++) ff2.push_back(ff[(i-j+7)%7]);
									hh.push_back(ff2);
								}
								sort(hh.begin(),hh.end());
								S.insert(hh[0]);
							}
	cout<<n<<" "<<S.size()<<"\n";
	for(vector<int> f:S){
		for(int i=0;i<7;i++) cout<<f[i]<<" ";
		cout<<"\n";
	}
} 

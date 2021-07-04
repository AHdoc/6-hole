#include<iostream>
#include<cstdio>
#include<cmath>
using namespace std;

const int MAXN=30;

struct Tpoint{
	double x,y;
	Tpoint(){}
	Tpoint(double _x,double _y){x=_x; y=_y;}
	Tpoint operator -(const Tpoint &b)const{return Tpoint(x-b.x,y-b.y);}
	double operator ^(const Tpoint &b)const{return x*b.y-y*b.x;}
	double operator *(const Tpoint &b)const{return x*b.x+y*b.y;}
};
double norm(Tpoint p){return sqrt(p.x*p.x+p.y*p.y);}

int n;
Tpoint p[MAXN+1];

int main(){
	freopen("3477710.txt","w",stdout);
	cin>>n;
	for(int i=1;i<=n;i++){
		double xx,yy;
		cin>>xx>>yy;
		p[i]=Tpoint(xx,yy);
	}
	cout<<n<<"\n";
	int tot=0;
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				int x=++tot;
				if(((p[j]-p[i])^(p[k]-p[i]))>0) cout<<tot<<" ";
				else cout<<"-"<<tot<<" ";
			}
	cout<<"0\n";
}

/*
[8730]
18
125 171
142 209
169 172
109 227
177 224
94 151
149 135
0 208
0 165
201 206
196 279
90 268
312 151
212 0
88 10
305 239
200 161
87 198

[347771]
29
1 1260
16 743
22 531
37 0
306 592
310 531
366 552
371 487
374 525
392 575
396 613
410 539
416 550
426 526
434 552
436 535
446 565
449 518
450 498
453 542
458 526
489 537
492 502
496 579
516 467
552 502
754 697
777 194
1259 320
*/

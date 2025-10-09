#include<iostream> 
#include<math.h> 
using namespace std;
#define R 8.31451 
double get_a(double w,double T,double Tc,double pc) 
{ 
double k=0.37464+1.54226*w-0.26992*w*w; 
double ar=(1+k*(1-sqrt(T/Tc)))*(1+k*(1-sqrt(T/Tc))); 
double a=0.45724*ar*R*R*Tc*Tc/pc; 
return a; 
} 
double get_b(double Tc,double pc) 
{ 
double b=0.0778*R*Tc/pc; 
return b; 
}
double Newton(double A,double B,double x) 
{ 
double x0; 
double f,df; 
do 
{ 
x0=x; 
f=x*x*x-(1-B)*x*x+(A-3*B*B-2*B)*x-(A*B-B*B-B*B*B); 
df=3*x*x-2*(1-B)*x+(A-3*B*B-2*B); 
x=x-f/df; 
}while(fabs(x-x0)>1e-6); 
return x; 
} 
double get_ar(double T,double v,double vv,double a,double b) 
{ 
double ar=R*T*log((v-b)/v)-a/(2*1.414*b)*(log((v-0.414*b)/(v+2.414*b))) +R*T*log(v/vv); 
return ar; 
} 
double get_sr(double T,double v,double vv,double a,double b,double bb) 
{ 
double sr=-1*R*log((v-b)/v)+bb/(2*1.414*b)*(log((v-0.414*b)/(v+2.414*b))) -R*log(v/vv); 
return sr; 
} 
void get_hr(double Tc,double pc,double w,double T,double p,double *hr,double *sr,double zz) 
{ 
double a=get_a(w,T,Tc,pc); 
double b=get_b(Tc,pc); 
double bb=(get_a(w,T+0.25,Tc,pc)-get_a(w,T-0.25,Tc,pc))/0.5; 
double A=a*p/(R*R*T*T); 
double B=b*p/(R*T); 
double z=Newton(A,B,zz); 
double v=z*R*T/p; 
double vv=R*T/p; 
double ar=get_ar(T,v,vv,a,b); 
*sr=get_sr(T,v,vv,a,b,bb); 
*hr=ar+T*(*sr)+R*T*(1-z); 
} 
int main() 
{ 
int i; 
double M; 
double w;
double h0=200.0; 
double s0=1.0; 
double T0; 
double p0; 
double c0,c1,c2,c3; 
double T,p,Tc,pc; 
double hr0,sr0,hrv,hrl,srv,srl; 
double pM[4]={44.096,58.122,114.04,114.04}; 
double pTc[4]={369.89,407.81,367.85,382.52}; 
double ppc[4]={4251200,3629000,3382200,3636300}; 
double pT0[4]={273.15,273.15,273.15,273.15}; 
double pp0[4]={474460,156960,315820,216480}; 
double pw[4]={0.1521,0.184,0.276,0.313}; 
double pc0[4]={-95.80,-23.91,18.349,55.389}; 
double pc1[4]={6.945,6.605,128.316,10.784}; 
double pc2[4]={-3.597*1e-3,-3.176*1e-3,-33.354,99.250}; 
double pc3[4]={7.290*1e-7,-4.981*1e-7,2.086,-49.88};
N1: 
cout<<"please enter 1(R290),2(R600a),3(R1234yf),4(R1234ze(E))"<<endl; 
cin>>i; 
if(i==1||i==2||i==3||i==4) 
{ 
M=pM[i-1]; 
Tc=pTc[i-1]; 
pc=ppc[i-1]; 
T0=pT0[i-1]; 
p0=pp0[i-1]; 
w=pw[i-1]; 
c0=pc0[i-1]; 
c1=pc1[i-1]; 
c2=pc2[i-1]; 
c3=pc3[i-1]; 
} 
else 
{ 
cout<<"The number is wrong"<<endl; 
goto N1; 
} 
cout<<"please enter T(K)"<<endl; 
cin>>T; 
cout<<"please enter p(Mpa)"<<endl; 
cin>>p;
p=p*1e6; 
get_hr(Tc,pc,w,T0,p0,&hr0,&sr0,0.001); 
get_hr(Tc,pc,w,T,p,&hrv,&srv,1.1); 
get_hr(Tc,pc,w,T,p,&hrl,&srl,0.001); 
if(fabs(hrv-hrl)<1e-4) 
{ 
double h=h0*M+hr0-hrv+R*(c0*(T-T0)+c1/2/Tc*(pow(T,2) 
-pow(T0,2))+c2/3/pow(Tc,2)*(pow(T,3)-pow(T0,3))+c3/4/pow(Tc,3) 
*(pow(T,4)-pow(T0,4))); 
double s=s0*M+sr0-srv-R*log(p/p0)+R*(c0*log(T/T0)+c1*(T-T0)/Tc 
+c2/2/pow(Tc,2)*(pow(T,2)-pow(T0,2))+c3/3/pow(Tc,3)*(pow(T,3)-pow(T0,3))); 
h=h/M; 
s=s/M; 
cout<<"h="<<h<<"kJ/kg"<<endl; 
cout<<"s="<<s<<"kJ/(kg*K)"<<endl; 
} 
else 
{ 
double h=h0*M+hr0-hrv+R*(c0*(T-T0)+c1/2/Tc*(pow(T,2) 
-pow(T0,2))+c2/3/pow(Tc,2)*(pow(T,3)-pow(T0,3))+c3/4/pow(Tc,3) 
*(pow(T,4)-pow(T0,4))); 
double s=s0*M+sr0-srv-R*log(p/p0)+R*(c0*log(T/T0)+c1*(T-T0)/Tc 
+c2/2/pow(Tc,2)*(pow(T,2)-pow(T0,2))+c3/3/pow(Tc,3)*(pow(T,3)-pow(T0,3))); 
h=h/M; 
s=s/M; 
cout<<" ÆøÏà "<<endl; 
cout<<"h="<<h<<"kJ/kg"<<endl; 
cout<<"s="<<s<<"kJ/(kg*K)"<<endl; 
h=h0*M+hr0-hrl+R*(c0*(T-T0)+c1/2/Tc*(pow(T,2) 
-pow(T0,2))+c2/3/pow(Tc,2)*(pow(T,3)-pow(T0,3))+c3/4/pow(Tc,3) 
*(pow(T,4)-pow(T0,4))); 
s=s0*M+sr0-srl-R*log(p/p0)+R*(c0*log(T/T0)+c1*(T-T0)/Tc 
+c2/2/pow(Tc,2)*(pow(T,2)-pow(T0,2))+c3/3/pow(Tc,3)*(pow(T,3)-pow(T0,3))); 
h=h/M; 
s=s/M; 
cout<<" ÒºÏà "<<endl; 
cout<<"h="<<h<<"kJ/kg"<<endl; 
cout<<"s="<<s<<"kJ/(kg*K)"<<endl; 
} 
}



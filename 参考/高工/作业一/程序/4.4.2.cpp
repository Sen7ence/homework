#include<iostream> 
#include<math.h> 
using namespace std;
#define R 8.31451 
#define k12 0.01 
double get_a(double w1,double w2,double T,double Tc1,double pc1,double Tc2, 
double pc2,double x1,double x2) 
{ 
double k=0.37464+1.54226*w1-0.26992*w1*w1; 
double ar=(1+k*(1-sqrt(T/Tc1)))*(1+k*(1-sqrt(T/Tc1))); 
double a1=0.45724*ar*R*R*Tc1*Tc1/pc1; 
k=0.37464+1.54226*w2-0.26992*w2*w2; 
ar=(1+k*(1-sqrt(T/Tc2)))*(1+k*(1-sqrt(T/Tc2))); 
double a2=0.45724*ar*R*R*Tc2*Tc2/pc2; 
double a=2*x1*x2*(1-k12)*sqrt(a1*a2)+x1*x1*a1+x2*x2*a2; 
return a; 
} 
double get_b(double Tc1,double pc1,double Tc2,double pc2,double x1,double x2) 
{ 
double b1=0.0778*R*Tc1/pc1; 
double b2=0.0778*R*Tc2/pc2; 
double b=x1*b1+x2*b2; 
return b; 
} 
double get_bb(double w1,double w2,double T,double Tc1,double pc1,double Tc2, 
double pc2,double x1,double x2) 
{ 
double a1=get_a(w1,w2,T+0.25,Tc1,pc1,Tc2,pc2,x1,x2); 
double a2=get_a(w1,w2,T-0.25,Tc1,pc1,Tc2,pc2,x1,x2); 
double bb=(a1-a2)/0.5; 
return bb; 
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
double ar=R*T*log((v-b)/v)-a*log((v-0.414*b)/(v+2.414*b))/(2*1.414*b)+R*T*log(v/vv); 
return ar; 
} 
double get_sr(double T,double v,double vv,double a,double b,double bb) 
{ 
double sr=-1*R*log((v-b)/v)+bb*log((v-0.414*b)/(v+2.414*b))/(2*1.414*b)-R*log(v/vv); 
return sr; 
} 
int main() 
{  
double M1=44.096; 
double Tc1=369.89; 
double pc1=4251200; 
double Ts1=273.15; 
double ps1=474460; 
double w1=0.1521; 
double x1; 
double c01=-95.80; 
double c11=6.945; 
double c21=-3.597*1e-3; 
double c31=7.290*1e-7; 
double M2=58.122; 
double Tc2=407.81; 
double pc2=3629000; 
double Ts2=273.15; 
double ps2=156960; 
double w2=0.184; 
double x2; 
double c02=-23.91; 
double c12=6.605; 
double c22=-3.176*1e-3; 
double c32=4.981*1e-7; 
double T0=273.15; 
double p0=329790; 
double hr0,sr0,hrv,srv,hrl,srl,h,s; 
double p,T; 
x1=(M2)/(M1+M2); 
x2=1-x1; 
double M=M1*x1+M2*x2; 
double a=get_a(w1,w2,T0,Tc1,pc1,Tc2,pc2,x1,x2); 
double b=get_b(Tc1,pc1,Tc2,pc2,x1,x2); 
double bb=get_bb(w1,w2,T0,Tc1,pc1,Tc2,pc2,x1,x2); 
double A=a*p0/(R*R*T0*T0); 
double B=b*p0/(R*T0); 
double z=Newton(A,B,0.001); 
double v=z*R*T0/p0; 
double vv=R*T0/p0; 
double ar=get_ar(T0,v,vv,a,b); 
sr0=get_sr(T0,v,vv,a,b,bb); 
hr0=ar+T0*sr0+R*T0*(1-z); 
cout<<"please enter T(K)"<<endl; 
cin>>T; 
cout<<"please enter p(Mpa)"<<endl; 
cin>>p; 
p=p*1e6; 
a=get_a(w1,w2,T,Tc1,pc1,Tc2,pc2,x1,x2); 
bb=get_bb(w1,w2,T,Tc1,pc1,Tc2,pc2,x1,x2); 
A=a*p/(R*R*T*T); 
B=b*p/(R*T); 
z=Newton(A,B,0.001);  
v=z*R*T/p; 
vv=R*T/p; 
ar=get_ar(T,v,vv,a,b); 
srv=get_sr(T,v,vv,a,b,bb); 
hrv=ar+T*srv+R*T*(1-z); 
z=Newton(A,B,1.1); 
v=z*R*T/p; 
vv=R*T/p; 
ar=get_ar(T,v,vv,a,b); 
srl=get_sr(T,v,vv,a,b,bb); 
hrl=ar+T*srl+R*T*(1-z); 
double dh1=R*(c01*(T-T0)+c11/2/Tc1*(pow(T,2)-pow(T0,2))+ 
c21/3/pow(Tc1,2)*(pow(T,3)-pow(T0,3))+c31/4/pow(Tc1,3)*(pow(T,4)-pow(T0,4))); 
double dh2=R*(c02*(T-T0)+c12/2/Tc2*(pow(T,2)-pow(T0,2))+ 
c22/3/pow(Tc2,2)*(pow(T,3)-pow(T0,3))+c32/4/pow(Tc2,3)*(pow(T,4)-pow(T0,4))); 
double dh=dh1*x1+dh2*x2; 
double ds1=-R*log(p/p0)+R*(c01*log(T/T0)+c11*(T-T0)/Tc1+ 
c21/2/pow(Tc1,2)*(pow(T,2)-pow(T0,2))+c31/3/pow(Tc1,3)*(pow(T,3)-pow(T0,3))); 
double ds2=-R*log(p/p0)+R*(c02*log(T/T0)+c12*(T-T0)/Tc2+ 
c22/2/pow(Tc2,2)*(pow(T,2)-pow(T0,2))+c32/3/pow(Tc2,3)*(pow(T,3)-pow(T0,3))); 
double ds=ds1*x1+ds2*x2; 
if(fabs(hrv-hrl)<1e-4) 
{ 
h=200+(hr0-hrv+dh)/M; 
s=1.0+(sr0-srv+ds)/M; 
cout<<"h="<<h<<"kJ/kg"<<endl; 
cout<<"s="<<s<<"kJ/(kg*K)"<<endl; 
} 
else 
{ 
h=200+(hr0-hrv+dh)/M; 
s=1.0+(sr0-srv+ds)/M; cout<<" ÆøÏà "<<endl; 
cout<<"h="<<h<<"kJ/kg"<<endl; 
cout<<"s="<<s<<"kJ/(kg*K)"<<endl; 
h=200+(hr0-hrl+dh)/M; 
s=1.0+(sr0-srl+ds)/M; cout<<" ÒºÏà "<<endl; cout<<"h="<<h<<"kJ/kg"<<endl; 
cout<<"s="<<s<<"kJ/(kg*K)"<<endl; 
}  
} 

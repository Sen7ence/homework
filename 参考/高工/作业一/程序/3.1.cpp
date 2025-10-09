#include<iostream> 
#include "math.h" 
using namespace std;
#define R 8.31451 
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
void R290(double T,double p,double *a,double *b,double *M) 
{ 
double Tc,pc,w,k,a1,Tr; 
*M=44.096e-3; 
Tc=369.89; 
pc=4251200; 
w=0.1521; 
k=0.37464+1.54226*w-0.26992*w*w; 
Tr=T/Tc; 
a1=pow(1+k*(1-pow(Tr,0.5)),2); 
*a=0.45727*a1*R*R*Tc*Tc/pc; 
*b=0.07780*R*Tc/pc; 
} 
void R600a(double T,double p,double *a,double *b,double *M) 
{ 
double Tc,pc,w,k,a1,Tr;
*M=58.122e-3; 
Tc=407.81; 
pc=3629000; 
w=0.184; 
k=0.37464+1.54226*w-0.26992*w*w; 
Tr=T/Tc; 
a1=pow(1+k*(1-pow(Tr,0.5)),2); 
*a=0.45727*a1*R*R*Tc*Tc/pc; 
*b=0.07780*R*Tc/pc; 
} 
void R1234yf(double T,double p,double *a,double *b,double *M) 
{ 
double Tc,pc,w,k,a1,Tr;
*M=114.04e-3; 
Tc=367.85; 
pc=3382200; 
w=0.276; 
k=0.37464+1.54226*w-0.26992*w*w; 
Tr=T/Tc; 
a1=pow(1+k*(1-pow(Tr,0.5)),2); 
*a=0.45727*a1*R*R*Tc*Tc/pc; 
*b=0.07780*R*Tc/pc; 
} 
void R1234ze(double T,double p,double *a,double *b,double *M) 
{ 
double Tc,pc,w,k,a1,Tr;
*M=114.04e-3; 
Tc=382.52; 
pc=3636300; 
w=0.313; 
k=0.37464+1.54226*w-0.26992*w*w; 
Tr=T/Tc; 
a1=pow(1+k*(1-pow(Tr,0.5)),2); 
*a=0.45727*a1*R*R*Tc*Tc/pc; 
*b=0.07780*R*Tc/pc; 
} 
void Hun(double T,double p,double *a,double *b,double *M) 
{ 
double a1,a2,b1,b2,x1,x2,k12,M1,M2; 
k12=0.01; 
R290(T,p,&a1,&b1,&M1); 
R600a(T,p,&a2,&b2,&M2); 
x1=1/(1+M1/M2); 
x2=1/(1+M2/M1); 
*a=x1*x1*a1+x2*x2*a2+2*x1*x2*(1-k12)*sqrt(a1*a2); 
*b=x1*b1+x2*b2; 
*M=x1*M1+x2*M2; 
} 
int main() 
{ 
double M,T,a,b,p,A,B; 
int i; 
N1: 
cout<<"please enter 1(R290),2(R600a),3(R1234yf),4(R1234ze)or5(Hun)"<<endl; 
cin>>i; 
if(i!=1&&i!=2&&i!=3&&i!=4&&i!=5) 
{ 
cout<<"The number is wrong"<<endl; 
goto N1; 
} 
cout<<"please enter T(K)"<<endl; 
cin>>T; 
cout<<"please enter p(Mpa)"<<endl; 
cin>>p; 
p=p*1e6; 
if(i==1) 
{ 
R290(T,p,&a,&b,&M); 
} 
else if(i==2)
{ 
R600a(T,p,&a,&b,&M); 
} 
else if(i==3) 
{ 
R1234yf(T,p,&a,&b,&M); 
} 
else if(i==4) 
{ 
R1234ze(T,p,&a,&b,&M); 
}
 else if(i==5) 
{ 
Hun(T,p,&a,&b,&M); 
} 
A=a*p/(R*R*T*T); 
B=b*p/(R*T); 
double z1=Newton(A,B,1000); 
double z2=Newton(A,B,0.001); 
if(fabs(z1-z2)<1e-4) 
{ 
double v1=z1*R*T/p/M; 
cout<<" 单位比体积为："<<v1<<"m^3/kg"<<endl; 
} 
else 
{ 
double v1=z1*R*T/p/M; 
double v2=z2*R*T/p/M; 
cout<<" 气体单位比体积为："<<v1<<"m^3/kg"<<endl; 
cout<<" 液体单位比体积为："<<v2<<"m^3/kg"<<endl; 
} 
} 

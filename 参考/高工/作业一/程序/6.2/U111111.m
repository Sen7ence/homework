Tc=[369.89,407.81];
Pc=[4.2512,3.629]*10^6;
w=[0.1521,0.184];
k12=0.01;
R=8.31451;%J/(mol·K)
k=0.37464+1.54226*w-0.26992*w.^2;
b=0.07780*R.*Tc./Pc;
fp1=zeros(1,2);fp2=zeros(1,2);
u= input('请输入压力序号: 1--10atm 2--1atm ');
for y1=0:0.001:1
if u==1
T=295;P=0.101e7; %P=10atm
else
T=213;P=0.101e6;%P=1atm
end
y2=1-y1;x1=0.1;x2=1-x1;x=0;
while abs(x-1)>=1e-2
T=T+0.1;
m=(1+k.*(1-(T./Tc).^0.5)).^2;
a=0.45724*m*R^2.*Tc.^2./Pc;
for i=1:2
if i==1
Z=0.0001;
am=x1^2*a(1,1)+2*x1*x2*(1-k12)*sqrt(a(1,1)*a(1,2))+x2^2*a(1,2);
bm=x1*b(1,1)+x2*b(1,2);
else
Z=1.1; 
am=y1^2*a(1,1)+2*y1*y2*(1-k12)*sqrt(a(1,1)*a(1,2))+y2^2*a(1,2);
bm=y1*b(1,1)+y2*b(1,2);
end
A=am*P/(R*T)^2;
B=bm*P/R/T;
f=Z^3-(1-B)*Z^2+(A-3*B^2-2*B)*Z-(A*B-B^2-B^3);
f1=3*Z^2-2*(1-B)*Z+(A-3*B^2-2*B);
Y=Z-f/f1;
while abs(Y-Z)>10^(-6)
Z=Y;
f=Z^3-(1-B)*Z^2+(A-3*B^2-2*B)*Z-(A*B-B^2-B^3);
f1=3*Z^2-2*(1-B)*Z+(A-3*B^2-2*B);
Y=Z-f/f1;
end
if i==2
fp1(1,i)=exp(b(1,1)/bm*(Y-1)-log(Y-B)-A/B/sqrt(8)*(2*(y1*a(1,1)+y2*(1-k12)*sqrt(a(1,1)*a(1,2)))/am-b(1,1)/bm)*log((Y+2.414*B)/(Y-0.414*B)));
fp2(1,i)=exp(b(1,2)/bm*(Y-1)-log(Y-B)-A/B/sqrt(8)*(2*(y2*a(1,2)+y1*(1-k12)*sqrt(a(1,1)*a(1,2)))/am-b(1,2)/bm)*log((Y+2.414*B)/(Y-0.414*B)));
else
fp1(1,i)=exp(b(1,1)/bm*(Y-1)-log(Y-B)-A/B/sqrt(8)*(2*(x1*a(1,1)+x2*(1-k12)*sqrt(a(1,1)*a(1,2)))/am-b(1,1)/bm)*log((Y+2.414*B)/(Y-0.414*B)));
fp2(1,i)=exp(b(1,2)/bm*(Y-1)-log(Y-B)-A/B/sqrt(8)*(2*(x2*a(1,2)+x1*(1-k12)*sqrt(a(1,1)*a(1,2)))/am-b(1,2)/bm)*log((Y+2.414*B)/(Y-0.414*B)));
end
end
k1=fp1(1,2)/fp1(1,1);
k2=fp2(1,2)/fp2(1,1);
x1=k1*y1/(k1*y1+k2*y2);
x2=k2*y2/(k1*y1+k2*y2);
x0=x;
x=k1*y1+k2*y2;
while abs(x-x0)>1e-3
Z=0.0001;
am=x1^2*a(1,1)+2*x1*x2*(1-k12)*sqrt(a(1,1)*a(1,2))+x2^2*a(1,2);
bm=x1*b(1,1)+x2*b(1,2);
A=am*P/(R*T)^2;
B=bm*P/R/T;
f=Z^3-(1-B)*Z^2+(A-3*B^2-2*B)*Z-(A*B-B^2-B^3);
f1=3*Z^2-2*(1-B)*Z+(A-3*B^2-2*B);
Y=Z-f/f1;
while abs(Y-Z)>10^(-6)
Z=Y;
f=Z^3-(1-B)*Z^2+(A-3*B^2-2*B)*Z-(A*B-B^2-B^3);
f1=3*Z^2-2*(1-B)*Z+(A-3*B^2-2*B);
Y=Z-f/f1;
end
fp1(1,1)=exp(b(1,1)/bm*(Y-1)-log(Y-B)-A/B/sqrt(8)*(2*(x1*a(1,1)+x2*(1-k12)*sqrt(a(1,1)*a(1,2)))/am-b(1,1)/bm)*log((Y+2.414*B)/(Y-0.414*B)));
fp2(1,1)=exp(b(1,2)/bm*(Y-1)-log(Y-B)-A/B/sqrt(8)*(2*(x2*a(1,2)+x1*(1-k12)*sqrt(a(1,1)*a(1,2)))/am-b(1,2)/bm)*log((Y+2.414*B)/(Y-0.414*B)));
k1=fp1(1,2)/fp1(1,1);
k2=fp2(1,2)/fp2(1,1);
x1=k1*y1/(k1*y1+k2*y2);
x2=k2*y2/(k1*y1+k2*y2);
x0=x;
x=k1*y1+k2*y2;
end
end
plot(x1,T,'b.')
hold on
plot(y1,T,'k.')
hold on
end
if u==1
title('溶液R290/R600a在10atm下的相图');
xlabel('组分x(y)');
ylabel('温度 T/K');
else
title('溶液R290/R600在1atm下的T-x相图');
xlabel('组分x(y)');
ylabel('温度 T/K');
end

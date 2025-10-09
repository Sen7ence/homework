p1=3e5;dp=100;
N=20000;err=1e-8;
u=1;
R=8.3145;
M1=44.096e-3;
Tc1=369.89;
Pc1=4.2512e6;
w1=0.1512;
k1=0.37464+1.54226*w1-0.26992*w1^2;
erre=1e-6;
for T1=200:0.1:369.89
for n=1:N
Tr1=T1/Tc1;
al1=(1+k1*(1-Tr1^0.5))^2;
a1=0.45724*al1*(R^2)*(Tc1^2)/Pc1;
b1=0.07780*R*Tc1/Pc1;
A=a1*p1/((R^2)*(T1^2));
B=b1*p1/(R*T1);
Z=1.1;
for o=0:1000
f=Z^3-(1-B)*Z^2+Z*(A-2*B-3*B^2)-(A*B-B^2-B^3);
Z=Z-f/(3*Z^2-2*(1-B)*Z+(A-2*B-3*B^2));
if(abs(f)<erre)
break
end
end
gongshi1v=exp(Z-1-log(Z-B)-A*log((Z+2.414*B)/(Z-0.414*B))/(2*sqrt(2)*B));
Z=0.001;
for p=0:1000
f=Z^3-(1-B)*Z^2+Z*(A-2*B-3*B^2)-(A*B-B^2-B^3);
Z=Z-f/(3*Z^2-2*(1-B)*Z+(A-2*B-3*B^2));
if(abs(f)<erre)
break
end
end
gongshi1L=exp(Z-1-log(Z-B)-A*log((Z+2.414*B)/(Z-0.414*B))/(2*sqrt(2)*B));
if abs(gongshi1v-gongshi1L)<=err
Y(u)=p1;
X(u)=T1;
u=u+1;
break
else
p1=p1+dp;
end
end
if n==N+1
fprintf('error!')
break;
else
hold on;
end
end
plot(X,Y/10^6,'r-');
grid;
title('R290a ¹¤ÖÊ p-T ÏàÍ¼');
xlabel('T/K');ylabel('p/MPa');

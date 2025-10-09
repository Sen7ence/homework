p1=3e5;p2=3e5;p3=3e5;p4=3e5;dp=100;
N=20000;er=1e-8;
l=1;
i=input('选择工质编号： 1(R290) 2(R600a) 3(R1234yf) 4(R1234ze(E))');
while i<5
switch i
case 1
for T1=200:0.01:369.89
for n=1:N
func1v=func1(T1,p1,1.1);
func1L=func1(T1,p1,0.001);
if abs(func1v-func1L)<=er
Y(l)=p1;
X(l)=T1;
l=l+1;
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
title('R290a工质的p-T相图');
xlabel('T/K');ylabel('p/MPa');
case 2
l=1;
for T2=200:0.01:369.89
for n=1:N
func2a=func2(T2,p2,1.1);
func2b=func2(T2,p2,0.001);
if abs(func2a-func2b)<=er
Y(l)=p2;
X(l)=T2;
l=l+1;
break
else
p2=p2+dp;
end
end
if n==N+1
fprintf('error!')
break;
else
hold on;
end
end
plot(X,Y/10^6,'b-');
grid;
title('R600a工质的p-T相图');
xlabel('T/K');ylabel('p/MPa');
case 3
l=1;
for T3=200:0.01:369.89
for n=1:N
func3a=func3(T3,p3,1.1);
func3b=func3(T3,p3,0.001);
if abs(func3a-func3b)<=er
Y(l)=p3;
X(l)=T3;
l=l+1;
break
else
p3=p3+dp;
end
end
if n==N+1
fprintf('error!')
break;
else
hold on;
end
end
plot(X,Y/10^6,'b-');
grid;
title('R1234yf工质的p-T相图');
xlabel('T/K');ylabel('p/MPa');
case 4
l=1;
for T4=200:0.01:369.89
for n=1:N
func4a=func4(T4,p4,1.1);
func4b=func4(T4,p4,0.001);
if abs(func4a-func4b)<=er
Y(l)=p4;
X(l)=T4;
l=l+1;
break
else
p4=p4+dp;
end
end
if n==N+1
fprintf('error!')
break;
else
hold on;
end
end
plot(X,Y/10^6,'k-');
grid;
title('R1234ze(E)工质的p-T相图');
xlabel('T/K');ylabel('p/MPa');
end
i=input('选择工质编号: 1--R290 2--R600a 3--R1234yf 4--R1234ze(E)');
end

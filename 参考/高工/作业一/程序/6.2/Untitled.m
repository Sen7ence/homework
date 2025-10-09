Tc=647.14;
Pc=22.046e6;
e=2.71828;
for i=1:547.14  %¸
    T(i)=100+1*i;
    Tr(i)=T(i)/Tc;
    Pr(i)=exp(log(Tr(i))*(7.60794067+10.1932439*(1-Tr(i))^1.89+21.1083545*(1-Tr(i))^5.67));
    P(i)=Pr(i)*Pc;
end   
plot(T,P,'m -'); 
grid;
title('水相变p-T图');
xlabel('T/K')
ylabel('P/Pa')
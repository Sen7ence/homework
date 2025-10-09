function Z=new(A,B,Z)
er=1e-6;
for n=0:1000
f=Z^3-(1-B)*Z^2+Z*(A-2*B-3*B^2)-(A*B-B^2-B^3);
Z=Z-f/(3*Z^2-2*(1-B)*Z+(A-2*B-3*B^2));
if(abs(f)<er)
break
end
end
end
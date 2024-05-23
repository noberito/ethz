function [s1,s2]=Principal(s)
m=(s(1)+s(2))/2;
d=sqrt((s(1)-s(2))^2+4*s(3)^2)/2;
s1=m+d;
s2=m-d;
end
function peeq = EMC(eta,theta,a,b,c,n)
global eps sig;
%
%Hosford criterion
f1 =  2/3*cos(pi/6*(1-theta));
f2 =  2/3*cos(pi/6*(3+theta));
f3 = -2/3*cos(pi/6*(1+theta));
Hf=(0.5*(f1-f3).^a+0.5*(f1-f2).^a+0.5*(f2-f3).^a).^(1/a);
peeq =b.*(1+c)^(1/n)*(Hf+c.*(2*eta+f1+f3))^(-1/n);
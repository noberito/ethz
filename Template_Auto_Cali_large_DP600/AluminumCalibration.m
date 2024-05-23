function [alpha,P_f,error_f]=AluminumCalibration(Pi,a0)
global a0 P Pi Sigma K0;

P=zeros(length(Pi),1);
K0=Pi(1);
a=[1 1 1 1 1 1 1 1];
    
Sigma=[1/2;1/2;1/2];
[Phi,dPhi]=EvalPhi(100,a);
options = optimset('Display','off','TolFun',1e-15, 'TolX',1e-15,'MaxFunEvals',100000,'MaxIter',100000);
[A,fval,exitflag,output] = fminsearch(@CalculateError,a,options);

alpha=A;
P_f=[P Pi P./Pi]
error=fval;
error_f=error;
end
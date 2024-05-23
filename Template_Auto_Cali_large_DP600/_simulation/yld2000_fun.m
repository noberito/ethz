function [yld2000_val] = yld2000_fun(sigma, a, L1, L2)
%Evaluates the Value of the yld-2000-2d yield function for a given set of
%alpha and a given a at a given stress state sigma
%
%Input: 
%
%sigma = [sigma_xx; sigma_yy; sigma_xy]
%a     = 8 for fcc and = 6 for bcc
%
%Output:
%
%yld2000_val = Phi_alpha(sigma)



X1 = L1 * sigma;
X2 = L2 * sigma;

X1_prin1 = 1/2*( X1(1) + X1(2) + sqrt( (X1(1)-X1(2))^2 + 4*X1(3)^2));
X1_prin2 = 1/2*( X1(1) + X1(2) - sqrt( (X1(1)-X1(2))^2 + 4*X1(3)^2));

X2_prin1 = 1/2*( X2(1) + X2(2) + sqrt( (X2(1)-X2(2))^2 + 4*X2(3)^2));
X2_prin2 = 1/2*( X2(1) + X2(2) - sqrt( (X2(1)-X2(2))^2 + 4*X2(3)^2));

Phi1 = abs(X1_prin1-X1_prin2)^a;
Phi2 = abs(2*X2_prin2+X2_prin1)^a + abs(2*X2_prin1+X2_prin2)^a;

yld2000_val = Phi1 + Phi2 - 2;

end

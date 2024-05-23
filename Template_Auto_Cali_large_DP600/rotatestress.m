global a0 Sigma;
a     = [0.7606;1.1813;1.0223;1.0028;1.0239;0.9224;1.0670;1.2060];
a0    = 8.0;

Sigma = [0;1;0];
for i = 1:length(S00(:,2))
[Phi,dPhidS] = EvalPhi(S00(i,2),a);
S90(i,1) = S00(i,2)^2/Phi;
end

Sigma = [1/2;1/2;1/2];
for i = 1:length(S00(:,2))
[Phi,dPhidS] = EvalPhi(S00(i,2),a);
S45(i,1) =  S00(i,2)^2/Phi;
end

plot(S00(:,1),S00(:,2),'b'); hold on
plot(S00(:,1),S90,'r'); hold on
plot(S00(:,1),S45,'g')
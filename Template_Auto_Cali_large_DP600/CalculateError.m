function error=CalculateError(a)

    global a0 Pi P Sigma K0;
    
    % Fit Parameters: Pi
    % Pi=[ Y0 r0 Y45 r45 Y90 r90 Ys0 Ys45];
  
    UT_Sigma=[1 0; 0 0];
    Mat_Ori=0:pi/2/2:pi/2;

for i=1:length(Mat_Ori)
   Q=[cos(Mat_Ori(i)) sin(Mat_Ori(i)); -sin(Mat_Ori(i)) cos(Mat_Ori(i))];
   Mat_Sigma=Q'*UT_Sigma*Q;
   Sigma=[Mat_Sigma(1,1);Mat_Sigma(2,2);Mat_Sigma(1,2)];
   [Phi,dPhi]=EvalPhi(100,a);
   Y(i)=100/Phi;
   Mat_PE=[dPhi(1) dPhi(3); dPhi(3) dPhi(2)];
   PE=Q*Mat_PE*Q';
   r(i)=-PE(2,2)/(PE(1,1)+PE(2,2));
end
for i=1:3
    P(i*2-1)=Y(i);
    P(2*i)=r(i);  
end

% Equi-biaxial°
% %      
Sigma=[1;1;0];
[Phi,dPhi]=EvalPhi(100,a);
YsEB=100/Phi;
rEB = dPhi(2)/dPhi(1);
P(7)=YsEB;
P(8)=rEB;  
       
% Uniaxial Shear, 0°
%     
% Sigma=[0;0;1];
% l=100;
% [Phi,dPhi]=EvalPhi(l,a);
% Ys0=100/Phi;
% P(15)=Ys0;
%     
% % Uniaxial Shear, 45°
%     
% Sigma=[1;-1;0];
% l=100;
% [Phi,dPhi]=EvalPhi(l,a);
% Ys45=l*K0/Phi;
% P(16)=Ys45;
    
   
    % Calculate the error term
%     for i=1:14
%         coeff=1.;
%         error(i)=((P(i)./Pi(i))-1.).^2*coeff;
%     end
    
    % Weight on Sig 0,45,90  
    for i=1:1:6
        coeff=10.; %10
        error(i)=((P(i)./Pi(i))-1).^2*coeff;
    end
%      % Weight on R 0,45,90  
%     for i=2:6:14
%         coeff=10.;
%         error(i)=((P(i)./Pi(i))-1).^2*coeff;
%     end
%     % Weight on rest Sig
%     for i=[3,5,9,11]
%         coeff=1.;
%         error(i)=((P(i)./Pi(i))-1).^2*coeff;
%     end
%     % Weight on rest R
%     for i=[4,6,10,12]
%         coeff=1.;
%         error(i)=((P(i)./Pi(i))-1).^2*coeff;
%     end
%     % Weight on rest EBT
    try
    for i=7:8
        coeff=2.;
        error(i)=((P(i)./Pi(i))-1).^2*coeff;
    end
    end
%     disp(['Error: ' num2str(error)]);
    error=sum(error);
    
end

function [Phi,dPhidS]=EvalPhi(lambda,a)

global a0 Sigma;

S=Sigma*lambda;

S1=[];
S1(1)=1/3*(2*a(1)*S(1)-a(1)*S(2));
S1(2)=1/3*(-a(2)*S(1)+2*a(2)*S(2));
S1(3)=a(7)*S(3);

S2=[];
S2(1)=1/9*((8*a(5)+2*(a(4)-a(6)-a(3)))*S(1)+(a(3)+4*(a(6)-a(4)-a(5)))*S(2));
S2(2)=1/9*((a(6)+4*(a(3)-a(4)-a(5)))*S(1)+(8*a(4)+2*(a(5)-a(6)-a(3)))*S(2));
S2(3)=a(8)*S(3);

[S1i,S1ii]=Principal(S1);
[S2i,S2ii]=Principal(S2);

Phi=1/(2^(1/a0))*((S1i-S1ii)^a0 + abs(2*S2i+S2ii)^a0 + abs(S2i+2*S2ii)^a0)^(1/a0);

% Calculation of derivatives, assuming a0=8 (a=8 in paper), 
% but should work for any even positive even integer.
dS1dS=zeros(3,3);
dS1dS(1,1)=2*a(1)/3;
dS1dS(1,2)=-a(1)/3;
dS1dS(1,3)=0;
dS1dS(2,1)=-a(2)/3;
dS1dS(2,2)=2*a(2)/3;
dS1dS(2,3)=0;
dS1dS(3,1)=0;
dS1dS(3,2)=0;
dS1dS(3,3)=a(7);


dS2dS=zeros(3,3);
dS2dS(1,1)=1/9*(8*a(5)+2*(a(4)-a(3)-a(6)));
dS2dS(1,2)=1/9*(a(3)+4*a(6)-4*a(4)-4*a(5));
dS2dS(1,3)=0;
dS2dS(2,1)=1/9*(a(6)+4*(a(3)-a(4)-a(5)));
dS2dS(2,2)=1/9*(8*a(4)+2*(a(5)-a(6)-a(3)));
dS2dS(2,3)=0;
dS2dS(3,1)=0;
dS2dS(3,2)=0;
dS2dS(3,3)=a(8);


dS1idS1=zeros(3,1);
Delta1=(S1(1)-S1(2))^2+4*S1(3)^2;
if Delta1>0
    dS1idS1(1)=1/2+(S1(1)-S1(2))/(2*sqrt(Delta1));
    dS1idS1(2)=1/2-(S1(1)-S1(2))/(2*sqrt(Delta1));
    % Note that S1_12 and S1_21 should be distinguished
    dS1idS1(3)=S1(3)/sqrt(Delta1);
end
dS1idS=zeros(3,1);
for i=1:3
    for j=1:3
        dS1idS(i)=dS1idS(i)+dS1idS1(j)*dS1dS(j,i);
    end
end

dS1iidS1=zeros(3,1);
Delta1=(S1(1)-S1(2))^2+4*S1(3)^2;
if Delta1>0
    dS1iidS1(1)=1/2-(S1(1)-S1(2))/(2*sqrt(Delta1));
    dS1iidS1(2)=1/2+(S1(1)-S1(2))/(2*sqrt(Delta1));
    % Note that S1_12 and S1_21 should be distinguished
    dS1iidS1(3)=-S1(3)/sqrt(Delta1);
end
dS1iidS=zeros(1,3);
for i=1:3
    for j=1:3
        dS1iidS(i)=dS1iidS(i)+dS1iidS1(j)*dS1dS(j,i);
    end
end


dS2idS2=zeros(1,3);
Delta2=(S2(1)-S2(2))^2+4*S2(3)^2;
if Delta2>0
    dS2idS2(1)=1/2+(S2(1)-S2(2))/(2*sqrt(Delta2));
    dS2idS2(2)=1/2-(S2(1)-S2(2))/(2*sqrt(Delta2));
    dS2idS2(3)=S2(3)/sqrt(Delta2);
end
dS2idS=zeros(1,3);
for i=1:3
    for j=1:3
        dS2idS(i)=dS2idS(i)+dS2idS2(j)*dS2dS(j,i);
    end
end

dS2iidS2=zeros(1,3);
Delta2=(S2(1)-S2(2))^2+4*S2(3)^2;
if Delta2>0
    dS2iidS2(1)=1/2-(S2(1)-S2(2))/(2*sqrt(Delta2));
    dS2iidS2(2)=1/2+(S2(1)-S2(2))/(2*sqrt(Delta2));
    dS2iidS2(3)=-S2(3)/sqrt(Delta2);
end
dS2iidS=zeros(1,3);
for i=1:3
    for j=1:3
        dS2iidS(i)=dS2iidS(i)+dS2iidS2(j)*dS2dS(j,i);
    end
end

dPhidS1i=Phi^(1-a0)/2*abs(S1i-S1ii)^(a0-1);
dPhidS1ii=-Phi^(1-a0)/2*(S1i-S1ii)^(a0-1);
dPhidS2i=Phi^(1-a0)/2*(2*abs(2*S2i+S2ii)^(a0-1)*sign(2*S2i+S2ii) ...
                       + abs(S2i+2*S2ii)^(a0-1)*sign(S2i+2*S2ii));
dPhidS2ii=Phi^(1-a0)/2*((2*S2i+S2ii)^(a0-1)*sign(2*S2i+S2ii) ...
                    + 2*(S2i+2*S2ii)^(a0-1)*sign(S2i+2*S2ii));

dPhidS=zeros(3,1);
for i=1:3
    dPhidS(i)=dPhidS1i*dS1idS(i)+dPhidS1ii*dS1iidS(i)+dPhidS2i*dS2idS(i)+dPhidS2ii*dS2iidS(i);
end

% if Sigma(1)+Sigma(2)==2*Sigma(3)
% Value=[S(3) S1(3) S2(3) a(7) a(8)]
% Deriv=[dS1dS(3,3) dS2dS(3,3) dS1idS1(3) dS1iidS1(3) dS2idS2(3) dS2iidS2(3)]
% dPhidSi=[dPhidS1i dPhidS1ii dPhidS2i dPhidS2ii]
% dPhi=dPhidS
% end

end

function [s1,s2]=Principal(s)
m=(s(1)+s(2))/2;
d=sqrt((s(1)-s(2))^2+4*s(3)^2)/2;
s1=m+d;
s2=m-d;
end


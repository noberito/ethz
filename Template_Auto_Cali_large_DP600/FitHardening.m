function [Swift,Voce,E_Swift,E_Voce] = FitHardening(Strain,Stress)
%Fit the Swift and Voce curves for different parts of the strain-stress data 
% Starting point : 0%,  1%,  5%,  10% of the total number of points
% End point      :90%, 95%, 99%, 100% of the total number of points

Length = length(Strain);
Start  = [1, floor(0.05*Length), floor(0.1*Length),floor(0.2*Length)];
End    = [floor(0.90*Length), floor(0.95*Length), floor(0.99*Length), Length];
% Start  = [1, floor(0.02*Length), floor(0.05*Length), floor(0.1*Length),floor(0.2*Length)];
% End    = [floor(0.70*Length), floor(0.8*Length),floor(0.90*Length), floor(0.95*Length), floor(0.99*Length), Length];

%Initialize results table and error table
Voce    = zeros(length(Start),length(End),3);
Swift   = zeros(length(Start),length(End),3);
E_Swift = zeros(length(Start),length(End),2);
E_Voce  = zeros(length(Start),length(End),2);

for i=1:length(Start)
    for j=1:length(End)
        Voce(i,j,:)  = FitGeneral(Strain(Start(i):End(j)),Stress(Start(i):End(j)),'Voce');
        Swift(i,j,:) = FitGeneral(Strain(Start(i):End(j)),Stress(Start(i):End(j)),'Swift');
        
        Swift_Stress = Swift_hardening(Strain,Swift(i,j,:));
        Voce_Stress  = Voce_hardening(Strain,Voce(i,j,:));
        E_Swift(i,j,1) = 1/length(Swift_Stress)*norm(Swift_Stress./Stress-1,1); %Compute the error for the whole curve
        E_Voce(i,j,1)  = 1/length(Voce_Stress)*norm(Voce_Stress./Stress-1,1);   %Compute the error for the whole curve
        
        
        StartError = floor(0.2*Length); %Compute the error in the center of the curve
        EndError   = floor(0.8*Length); %Compute the error in the center of the curve
        E_Swift(i,j,2) = 1/length(Swift_Stress(StartError:EndError))*norm(Swift_Stress(StartError:EndError)./Stress(StartError:EndError)-1,1);
        E_Voce(i,j,2)  = 1/length(Voce_Stress(StartError:EndError))*norm(Voce_Stress(StartError:EndError)./Stress(StartError:EndError)-1,1);

    end
end

end


function FitRes=FitGeneral(Strain,Stress,hardening)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%y = A*(x+e0)^n; Swift
%%%%y = k+Q*(1-exp(-b*x)); Voce
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch hardening
    case 'Swift'   
%                 1,    2,   3     
%   Swift    = [  A,   e0,   n];
    x0       = [0.6,   5.,  2.];
    Scale    = [1e3, 0.01, 0.1];
    lb = [0.1*x0(1),   0.,  0.]; %Lower bound
    ub = [ 10*x0(1),  20.,  4.]; %Upper bound
   
    case 'Voce'       
%                         1,    2,   3     
%   Voce  = [             k,    Q,   b];
    x0    = [Stress(1)/1000,   4.,  2.];
    Scale = [           1e3,  100,  10];
    lb = 0.1*x0;
    ub =  10*x0;    
    
end   

    %Fminsearch
    options = optimset('MaxFunEvals',1000,'MaxIter',1000,'TolFun',5e-4, 'TolX',5e-4);
    [x,fval,exitflag,output] = fminsearchbnd(@(param) ErrorFit(param,Strain,Stress,Scale,hardening), x0,lb,ub, options); % Nelder_Mead direct search (simplex, no derivatives)
    FitRes=x.*Scale;

end

function error = ErrorFit(param,Strain,Stress,Scale,hardening)
    weight = 0.98; %Weight for the maximum force criteria
    param=param.*Scale;
       
    switch hardening
        case 'Swift'  
            FitStress = Swift_hardening(Strain,param);
            Strain_neck = Swift_necking(param);
        case 'Voce'
            FitStress = Voce_hardening(Strain,param);
            Strain_neck  = Voce_necking(param);
    end
    
    Strain_neck_error = abs(Strain_neck/Strain(end)-1);  
    Curve_error = 1/length(FitStress)*norm(FitStress./Stress-1,1); %Calculate error    
    error = weight*Curve_error + (1-weight)*Strain_neck_error;
end

function Swift_Stress=Swift_hardening(Strain,param)
    Swift_Stress = param(1)*(Strain+param(2)).^param(3);
end

function Swift_neck=Swift_necking(param)
    Swift_neck = param(3)-param(2);
end

function Voce_Stress=Voce_hardening(Strain,param)
    Voce_Stress = param(1) + param(2)*(1-exp(-param(3)*Strain));
end

function Voce_neck=Voce_necking(param)
    Voce_neck = -1/param(3)*log((param(1)+param(2))/(param(2)*(1+param(3))));
end
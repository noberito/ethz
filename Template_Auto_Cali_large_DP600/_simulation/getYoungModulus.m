function Young = getYoungModulus(Strain,Stress) 

    %Fix the minimum point to start the interpolation    
    minStress = 0.08*max(Stress);
    minStrain = find(Stress<=minStress);
    minStrainInd = minStrain(end);
    minStrain = Strain(minStrainInd);
    
    %Fix the max strain 
    maxStrain = find(Strain<=0.01);
    maxStrainInd = maxStrain(end);
    maxStrain = Strain(maxStrainInd);
    
    
    TotalNumber = maxStrainInd - minStrainInd; %Total number of points
    StepSize = 5;                              %Stepsize for fitting
    NumberFit = floor(TotalNumber/StepSize);   %Total number of fit 
    
    %Initialize tables for A: Young, B: offset, R2: R-squared
    A  = nan(1,NumberFit);
    B  = nan(1,NumberFit);
    R2 = nan(1,NumberFit);
    
    %Start doing interpolation on larger and larger increment
    for i=1:NumberFit
        [A(i),B(i),R2(i)] = LinearRegression(Strain(minStrainInd:minStrainInd+i*StepSize),Stress(minStrainInd:minStrainInd+i*StepSize));
    end
    
    [~,I]=max(R2);
    Young = A(I); 
    
% %-------------------------------------------------------
% % Create PLOTS
% %-------------------------------------------------------     
%     plot(A)
%     yyaxis right
%     plot(R2)
% 
%     figure
%     
%     plot(Strain,Stress)
%     hold on
%     plot(Strain(1:maxStrainInd),Young*Strain(1:maxStrainInd)+B(I))

end

function [A,B,R2] = LinearRegression(X,Y) 

    %Fit linear polynom
    res = polyfit(X,Y,1);
   
    %Compute R-squared
    yfit = polyval(res,X);
    yresid = Y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(Y)-1) * var(Y);
    R2 = 1 - SSresid/SStotal;
    
    A = res(1);
    B = res(2);
    

end

function offset = FUoffset(displ,force) 

    %Fix the minimum point to start the interpolation    
    minStress = 0.08*max(force);
    minStrain = find(force<=minStress);
    minStrainInd = minStrain(end);
    minStrain = displ(minStrainInd);
    
    %Fix the max strain 
    maxStrain = find(displ<=0.1);
    maxStrainInd = maxStrain(end);
    maxStrain = displ(maxStrainInd);
    
    
    TotalNumber = maxStrainInd - minStrainInd; %Total number of points
    StepSize = 5;                              %Stepsize for fitting
    NumberFit = floor(TotalNumber/StepSize);   %Total number of fit 
    
    %Initialize tables for A: Young, B: offset, R2: R-squared
    A  = nan(1,NumberFit);
    B  = nan(1,NumberFit);
    R2 = nan(1,NumberFit);
    
    %Start doing interpolation on larger and larger increment
    for i=1:NumberFit
        [A(i),B(i),R2(i)] = LinearRegression(displ(minStrainInd:minStrainInd+i*StepSize),force(minStrainInd:minStrainInd+i*StepSize));
    end
    
    [~,I]=max(R2);
    offset = -B(I)/A(I);
    
% %-------------------------------------------------------
% % Create PLOTS
% % -------------------------------------------------------     
%     plot(A)
%     yyaxis right
%     plot(R2)
% 
%     figure
%     
%     plot(displ,force)
%     hold on
%     plot(displ(1:maxStrainInd),A(I)*displ(1:maxStrainInd)+B(I))

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

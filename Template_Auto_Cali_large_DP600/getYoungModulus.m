function Young = getYoungModulus(Strain,Stress,plotname) 
    %Fix the minimum point to start the interpolation    
    
    aa = find(Stress==max(Stress));
    aa = max(aa);
    Stress = Stress(1:aa,:);
    Strain = Strain(1:aa,:); 
    
   minStress = 0.15*max(Stress);
    minStrain = find(Stress<=minStress);
    minStrainInd = minStrain(end);
    minStrain = Strain(minStrainInd);
    
    %Fix the max strain 
    maxStrain = find(Strain<=0.008);
    maxStrainInd = maxStrain(end);
    maxStrain = Strain(maxStrainInd);
    
    
    TotalNumber = maxStrainInd - minStrainInd; %Total number of points
    %StepSize = 20;                              %Stepsize for fitting
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
    
%-------------------------------------------------------
% Create PLOTS
%-------------------------------------------------------     
    plot(A)
    yyaxis right
    plot(R2)

    figure(9)
    
    plot(Strain,Stress)
    hold on
    plot(Strain(1:maxStrainInd),Young*Strain(1:maxStrainInd)+B(I))
    pos=randi([length(Strain)-5,length(Stress)],1);
    eval(['text(Strain(pos),Stress(pos),''\leftarrow ' plotname ' '')']);
    set(gca,'XMinorTick','on');
    set(gca,'YMinorTick','on');
    set(gca,'TickLength',[0.01 0.01]);
    set(gca,'FontSize',18);
    set(gca,'LineWidth',1);
    % set(gca,'FontSize',18);
    % set(gca,'FontSize',18);
    xlabel('Eng. Strain [-]','FontSize',18);
    ylabel('Eng. Stress [MPa]','FontSize',18);
    axis square;
    box on

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

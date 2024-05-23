function Sim_pre_proc
clear all; close all; clc; format long
%-----------------------------------
% PARAMETERS
%-----------------------------------
minPStrain=0.001;
ratiostrain=0.08;
SVminstrain = 0.002;
%YLD2000
a0 = 6.; %exponent
YEBT = 1.0; %Yield Stress Ratio EBT
rEBT = 1.0; %"lankford" Ratio EBT
%-----------------------------------
% Load results from files
%-----------------------------------
cd _results
listOf = dir('*.csv');
ii=1;
for i=1:numel(listOf)
if (strncmp(listOf(i).name,'UT',2)==1 && strncmp(listOf(i).name,'UT_EBT',6)~=1)
filenames{ii}=listOf(i).name;
ii=ii+1;
end
end
clear i ii listOf
%-----------------------------------
PlasticStrain_int_max_all = 1.0;
for i=1:numel(filenames)
name=char(filenames(i))
geometry=  name(1:2);
orientation= name(4:5);
testnumber = name(7);
plotname = [geometry '-' orientation '-' testnumber];
testname = [geometry '_' orientation '_' testnumber];    
nameori = [geometry '_' orientation];
testnumber = str2num(name(7));
   
    A=csvread(char(filenames(i)),1,0);
    eval(['Lankford(testnumber).' nameori '=A(1,13);']);   
    
    a=A(:,10);
    aa=any(a~=0,2);
    A=A(aa,:);
    clear a aa
    
    [kk, ia, ic]  = unique(A(:,9));
    A=A(ia,:);
    clear kk ia ic
    [kk, ia, ic]  = unique(A(:,10));
    A=A(ia,:);
    clear kk ia ic
    
    eval(['PlasticStrain(testnumber).' nameori '=A(:,9);']);
    eval(['PlasticStress(testnumber).' nameori '=A(:,10);']);
  
%     eval(['PlasticStrain_max(testnumber)=max(PlasticStrain(testnumber).' nameori ');']);
       
    eval(['PlasticStrain_int(testnumber).' nameori '=(minPStrain:0.0001:A(end,9))'';']);
    eval(['PlasticStress_int(testnumber).' nameori '=pchip(PlasticStrain(testnumber).' nameori ',PlasticStress(testnumber).' nameori ',PlasticStrain_int(testnumber).' nameori ');']);
%     eval(['PlasticWork_int(testnumber).' nameori '=PlasticStrain_int(testnumber).' nameori '.*PlasticStress_int(testnumber).' nameori ';']);
    eval(['PlasticWork_int(testnumber).' nameori '=cumtrapz(PlasticStrain_int(testnumber).' nameori ', PlasticStress_int(testnumber).' nameori ');']);
    
    eval(['PlasticStrain_int_max_all = min(PlasticStrain_int_max_all,max(PlasticStrain_int(testnumber).' nameori '));']); 
        
    clear A a aa pos name geometry orientation plotname testname testnumber nameori

end
cd ..
%
PlasticStrain_int_all=(minPStrain:0.0001:PlasticStrain_int_max_all)';
ratiostrain = min(ratiostrain, PlasticStrain_int_max_all);
[minDistance, ratiolocation] = min(abs(PlasticStrain_int_all-ratiostrain));
% yieldstrain(1) = min(yieldstrain, PlasticStrain_int_max_all);
% [minDistance, yieldlocation] = min(abs(PlasticStrain_int_all-yieldstrain));
SVminstrain = min(SVminstrain, PlasticStrain_int_max_all);
[minDistance, SVminlocation] = min(abs(PlasticStrain_int_all-SVminstrain));
%
names=fieldnames(PlasticStress_int);
% SPLINE ALL SIG?EPSILON UP TO MAX EPS OF ALL
for i=1:numel(names);
name=char(names(i));   
geometry= name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];
for testnumber=1:100
try
eval(['PlasticStress_int_all(testnumber).' nameori '=pchip(PlasticStrain(testnumber).' nameori ',PlasticStress(testnumber).' nameori ',PlasticStrain_int_all);']);
eval(['PlasticWork_int_all(testnumber).' nameori '=pchip(PlasticStrain_int(testnumber).' nameori ',PlasticWork_int(testnumber).' nameori ',PlasticStrain_int_all);']);
end
end
clear i A a aa pos name geometry orientation plotname testname testnumber nameori
end
%-----------------------------------
% IDENTIFY OUTLIERS & CALCULATE MEAN
%-----------------------------------
for i=1:numel(names);
name=char(names(i));  
geometry=  name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];
for ii=1:100
try
eval(['temp=PlasticStress_int(ii).' nameori ';']);  
% eval(['SigY(ii).' nameori '=PlasticStress_int(ii).' nameori '(yieldlocation);']);
if temp>0
ll=ii;
end
end
clear temp ii
end
eval(['PlasticStress_int_mean.' nameori '=mean(cat(ll,PlasticStress_int_all.' nameori '),ll);']);
eval(['PlasticWork_int_mean.' nameori '=mean(cat(ll,PlasticWork_int_all.' nameori '),ll);']);
% DETERMINE STANDARD DEVIATION
eval(['tempStresses=squeeze(cat(ll,PlasticStress_int_all.' nameori '));']);
eval(['stdd.' nameori '=mean(std(tempStresses,0,2));']);
eval(['tempLank=squeeze(cat(ll,Lankford.' nameori '));']);
% DETERMINE ERROR
for ii=1:ll  
    eval(['err(ii)= mean(abs(PlasticStress_int_all(ii).' nameori ' - PlasticStress_int_mean.' nameori '));']);    
    if err(ii) > 2.*mean(std(tempStresses,0,2))
       usevector(1,ii)=0; 
    else
       usevector(1,ii)=1;
    end         
end 
eval(['[kk best.' nameori ']=min(err);']);
%CHOOSE EXP FOR SIM
eval(['PlasticStrain_sim.' nameori '=PlasticStrain(best.' nameori ').' nameori ';'])
eval(['PlasticStress_sim.' nameori '=PlasticStress(best.' nameori ').' nameori ';'])
eval(['PlasticStress_sim_int_all.' nameori '=tempStresses(:,best.' nameori ');'])
%CREATE MEAN
eval(['PlasticStress_int_mean_std.' nameori '=mean(tempStresses(:,logical(usevector)),2);']);
eval(['PlasticWork_int_all_mean_std.' nameori '=PlasticWork_int_mean.' nameori ';']);

%-----------------------------------
%  DETERMINE RATIO LOCATION ACC TO WORK EQUIVALENCE
%-----------------------------------
tempwork = PlasticWork_int_mean.UT_00;
yieldwork = tempwork(ratiolocation);
for ii=1:100
try
eval(['[yieldworklocation iib] = find(PlasticWork_int(ii).' nameori '<=yieldwork);']);
yieldworklocation=yieldworklocation(end)
eval(['SigY(ii).' nameori '=PlasticStress_int(ii).' nameori '(yieldworklocation);']);
end
clear temp ii yieldworklocation iib
end

eval(['[yieldworklocation iib] = find(PlasticWork_int_mean.' nameori '<=yieldwork);']);
yieldworklocation=yieldworklocation(end)
eval(['SigYRatio_mean.' nameori '=PlasticStress_int_mean_std.' nameori './PlasticStress_int_mean_std.UT_00;']);  
eval(['SigYRatio.' nameori '=SigYRatio_mean.' nameori '(yieldworklocation);']);
eval(['SigY_mean.' nameori '=mean([SigY(logical(usevector)).' nameori ']);']);  
eval(['Lank_mean.' nameori '=mean([Lankford(logical(usevector)).' nameori ']);']);

clear ll i ii name geometry orientation nameori tem* err
clear kk usevector yieldworklocation iib
end
disp('CHOSEN FOR SIMULATION:')
disp(best)
%-------------------------------------------------------
% Save data for Simulation
%------------------------------------------------------- 
mkdir _simulation
cd _simulation
for i=1:numel(names);
clear A    
name=char(names(i))   
geometry=  name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];

    eval(['A(:,1)=PlasticStrain_sim.' nameori ';']);
    eval(['A(:,2)=PlasticStress_sim.' nameori ';']);
    eval(['A(1,3)=best.' nameori ';']);
    
name = strcat(nameori,'_exp.csv')
dlmwrite(name,A);
clear A
end
cd ..
%-----------------------------------
% PLOT "AVERAGE TRUE STRESS-TRUE STRAIN"
%-----------------------------------
figure(3);
for i=1:numel(names);
name=char(names(i));   
geometry=  name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];
plotname =[geometry '-' orientation];
for ii=1:100
try
eval(['plot(PlasticStrain_int_all,PlasticStress_int_all(ii).' nameori ');']); hold on
pos=randi([length(PlasticStrain_int_all)-15,length(PlasticStrain_int_all)-5],1);
eval(['text(PlasticStrain_int_all(pos),PlasticStress_int_all(ii).' nameori '(pos),'' ' plotname ' '')']);
eval(['ymax(ii)=max(PlasticStress_int_all(ii).' nameori ');']); 
eval(['ymin(ii)=min(PlasticStress_int_all(ii).' nameori ');']);
end
end
eval(['plot(PlasticStrain_int_all,PlasticStress_int_mean_std.' nameori ',''.k'');']); hold on
eval(['plot(PlasticStrain_int_all,PlasticStress_sim_int_all.' nameori ',''.r'');']); hold on
% pos=randi([10,length(PlasticStrain_int_all)-5],1);
eval(['text(PlasticStrain_int_all(end),mean(PlasticStress_int_mean_std.' nameori '(end)),'' ' plotname ' '')']);
end
%SET IMAGE PROPERTIES
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
xmax=max(PlasticStrain_int_all);
ymin=min(ymin)-50;
ymax=max(ymax)+50;
eval(['axis([0 ' num2str(xmax) '  ' num2str(ymin) ' ' num2str(ymax) '])']);
xlabel('Plastic Strain [-]','FontSize',18);
ylabel('True Stress [MPa]','FontSize',18);
axis square;
box on
%-----------------------------------
% PLOT SIG RATIOS against plastic strain
%-----------------------------------
figure(1);
for i=1:numel(names);
name=char(names(i));   
geometry=  name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];
plotname =[geometry '-' orientation];
% eval(['SigYRatio_mean.' nameori '=PlasticStress_int_mean.' nameori './PlasticStress_int_mean.UT_00;']);
% eval(['SigYRatio.' nameori '=SigYRatio_mean.' nameori '(ratiolocation);']);
eval(['plot(PlasticStrain_int_all,SigYRatio_mean.' nameori ');']); hold on
pos=randi([10,length(PlasticStrain_int_all)-10],1);
eval(['text(PlasticStrain_int_all(pos),mean(PlasticStress_int_mean.' nameori '(pos-10:pos+10)./PlasticStress_int_mean.UT_00(pos-10:pos+10)),'' ' plotname ' '')']);
end
plot(PlasticStrain_int_all(ratiolocation),[0.5:0.002:1.5],'.r');
% plot(PlasticStrain_int_all(yieldlocation),[0.9:0.002:1.1],'.b');
plot(PlasticStrain_int_all(SVminlocation),[0.5:0.002:1.5],'.k');
%SET IMAGE PROPERTIES
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
xmax=max(PlasticStrain_int_all);
ymin=min(min(cell2mat(struct2cell(SigYRatio_mean))))-0.025;
ymax=max(max(cell2mat(struct2cell(SigYRatio_mean))))+0.025;
eval(['axis([0 ' num2str(xmax) '  ' num2str(ymin) ' ' num2str(ymax) '])']);
xlabel('Plastic Strain [-]','FontSize',18);
ylabel('\sigma_x/\sigma_0 [-]','FontSize',18);
axis square;
box on
%-----------------------------------
% PLOT SIG RATIOS against plastic WORK
%-----------------------------------
figure(10);
xmin = 50000;
xmax = 0;
for i=1:numel(names);
name=char(names(i));   
geometry=  name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];
plotname =[geometry '-' orientation];
eval(['PlasticWorkTemp = PlasticWork_int_all_mean_std.' nameori ';']);
% eval(['SigYRatio_mean.' nameori '=PlasticStress_int_mean.' nameori './PlasticStress_int_mean.UT_00;']);
% eval(['SigYRatio.' nameori '=SigYRatio_mean.' nameori '(ratiolocation);']);
% eval(['plot(PlasticStrain_int_all,SigYRatio_mean.' nameori ');']); hold on
eval(['plot(PlasticWorkTemp,SigYRatio_mean.' nameori ');']); hold on
pos=randi([10,length(PlasticWorkTemp)-10],1);
eval(['text(PlasticWorkTemp(pos),mean(PlasticStress_int_mean.' nameori '(pos-10:pos+10)./PlasticStress_int_mean.UT_00(pos-10:pos+10)),'' ' plotname ' '')']);
xmin=min(xmin,min(PlasticWorkTemp));
xmax=max(xmax,max(PlasticWorkTemp));
end
plot(PlasticWork_int_all_mean_std.UT_00(ratiolocation),[0.5:0.002:1.5],'.r');
% plot(PlasticStrain_int_all(yieldlocation),[0.9:0.002:1.1],'.b');
plot(PlasticWork_int_all_mean_std.UT_00(SVminlocation),[0.5:0.002:1.5],'.k');
%SET IMAGE PROPERTIES
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
ymin=min(min(cell2mat(struct2cell(SigYRatio_mean))))-0.025;
ymax=max(max(cell2mat(struct2cell(SigYRatio_mean))))+0.025;
eval(['axis([' num2str(xmin) ' ' num2str(xmax) '  ' num2str(ymin) ' ' num2str(ymax) '])']);
xlabel('Plastic Work [MJ/m^3]','FontSize',18);
ylabel('\sigma_x/\sigma_0 [-]','FontSize',18);
axis square;
box on
%-----------------------------------
% PLOT SIG & LANKFORD RATIOS against material orientation
%-----------------------------------
fig=figure(2);
set(fig,'color','w'); 
set(fig,'PaperSize', [5 5]);
hold on
x0=0;
y0=0;
width=600;
height=500;
set(fig,'units','points','position',[0,0,width,height])
%-----------------------------------
%Create first axis
%-----------------------------------
ymin=0.9*min(cell2mat(struct2cell(SigYRatio)));
ymax=1.1*max(cell2mat(struct2cell(SigYRatio)))+0.1;
eval(['axis([-5 95 ' num2str(ymin) ' ' num2str(ymax) '])']);
% axis([-5 95 0.8 1.2])
ax1=gca;
set(ax1,'Color','none');
set(ax1,'XMinorTick','on');
set(ax1,'XTick',[0:15:90]);
set(ax1,'YMinorTick','on');
set(ax1,'TickLength',[0.01 0.01]);
set(ax1,'FontSize',18);
set(ax1,'LineWidth',1);
xlabel('Material Orientation [°]','FontSize',18);
ylabel('\sigma_x_°/\sigma_0_° [-]','FontSize',18);
box(ax1,'off')
axis square
%-----------------------------------
%Create second axis
%-----------------------------------
ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','b','YColor','b');
linkaxes([ax1 ax2], 'x'); 
hold on
ymin=round(100*0.90*min(min(cell2mat(struct2cell(Lankford)))))/100-0.1;
ymax=round(100*1.1*max(max(cell2mat(struct2cell(Lankford)))))/100;
eval(['axis([-5 95 ' num2str(ymin) ' ' num2str(ymax) '])']);
% axis([-5 95 0.4 0.8])
axis square
set(ax2,'XTickLabel',[]);
set(ax2,'XMinorTick','on');
set(ax2,'XTick',[0:15:90],'fontsize',15);
% set(ax2,'YTick',[0:0.05:2.0],'fontsize',15);
set(ax2,'YMinorTick','on');
set(ax2,'TickLength',[0.01 0.01]);
set(ax2,'FontSize',18);
set(ax2,'LineWidth',1);
ylabel('r-Value [-]','fontsize',18,'color','b')
%-----------------------------------
%Plotting yield stress variation
%-----------------------------------
set(2,'CurrentAxes',ax1);
for i=1:numel(fieldnames(PlasticStress))
name=char(names(i))   
geometry=  name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];
plotname =[geometry '-' orientation];
eval(['h(1)=plot(str2num(orientation),SigYRatio.' nameori ',''.k'');']); hold on
%PLOT MAX/MIN
eval(['h(2)=plot(str2num(orientation),max(SigYRatio_mean.' nameori '),''sk'');']); hold on
eval(['h(3)=plot(str2num(orientation),min(SigYRatio_mean.' nameori '),''sk'');']); hold on
%PLOT ALL
%     eval(['plot(str2num(orientation),(SigYRatio_mean.' nameori '),''xk'');']); hold on
clear A a aa pos  
end
%-----------------------------------
% Plot R-Values & Save Plot
%-----------------------------------
set(2,'CurrentAxes',ax2);
for i=1:numel(fieldnames(Lankford))
name=char(names(i))   
geometry=  name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];
plotname =[geometry '-' orientation];    

    eval(['plot(str2num(orientation),Lank_mean.' nameori ',''.b'');']); hold on
    
for ii=1:100
try
eval(['temp=Lankford(ii).' nameori ';']);  
if temp>0
ll=ii;
end
end

end

%PLOT MAX/MIN
% eval(['plot(str2num(orientation),max(cat(ll,Lankford.' nameori ')),''xk'');']); hold on
% eval(['plot(str2num(orientation),min(cat(ll,Lankford.' nameori ')),''xk'');']); hold on
%PLOT ALL
eval(['plot(str2num(orientation),(cat(1,Lankford.' nameori ')),''xb'');']); hold on

clear ll ii name geometry orientation nameori temp
    
end
set(2,'CurrentAxes',ax1);
% axis square;
% h_l=legend([h(1) h(2)],{'Experiment','Yld2000-2d'},'location','SouthEast','FontSize',15);
cd _results
savefig(2,'__Y-R_ratios_Angle')
saveas(2,['__Y-R_ratios_Angle.png'],'png');
cd ..
clear plotname i xmax ymax ymin 
%-----------------------------------
% SAVE SigY and Lankford in Overview
%-----------------------------------
outputLank = struct2table(Lankford);
outputLankMean = struct2table(Lank_mean);
outputSigY = struct2table(SigY);
outputSigY_mean = struct2table(SigY_mean);
outputSigYRatio = struct2table(SigYRatio);

temp=table2array(struct2table(best));
outputSigTEpsPl = [PlasticStrain(temp(1)).UT_00, PlasticStress(temp(1)).UT_00];
colname = {'UT-0° Exp Plastic Strain','UT-0° Exp True Stress [MPa]'};
outputSigTEpsPl = array2table(outputSigTEpsPl,'VariableNames',colname);
clear temp

cd _results
fname = '_Experiments_Overview.xlsx'
writetable(outputSigY,fname,'Sheet','Yld Stress Ratio Stresses (MPa)');
writetable(outputLank,fname,'Sheet','Lankford Ratios');
writetable(outputSigTEpsPl,fname,'Sheet','Hardening');
% writetable(outputSigY_mean,fname,'Sheet','Mean Yield Stresses in MPa');
% writetable(outputSigYRatio,fname,'Sheet','Yield Stress Ratios');
% writetable(outputLankMean,fname,'Sheet','Mean Lankford Ratios');
clear fname output*
cd ..
%-----------------------------------
% SAVE PLOTS
%-----------------------------------
cd _results
savefig(1,'__Y_ratios_Strain')
saveas(1,['__Y_ratios_Strain.png'],'png');
savefig(10,'__Y_ratios_PlasticWork')
saveas(10,['__Y_ratios_PlasticWork.png'],'png');
savefig(2,'__Y-R_ratios_Angle')
saveas(2,['__Y-R_ratios_Angle.png'],'png');
savefig(3,'__Best_UTs')
saveas(3,['__Best_UTs.png'],'png');
cd ..
mkdir _simulation
cd _simulation
savefig(2,'__Y-R_ratios_Angle')
saveas(2,['__Y-R_ratios_Angle.png'],'png');
cd ..
%-------------------------------------------------------
% Fit SV Hardening
%-------------------------------------------------------    
aa=any(PlasticStrain_sim.UT_00>=SVminstrain,2);
Strain = PlasticStrain_sim.UT_00;
Strain = Strain(aa,:);
Stress = PlasticStress_sim.UT_00;
Stress = Stress(aa,:);

[S,V,E_S,E_V]=FitHardening(Strain,Stress);

minMatrix = min(min(E_S(:,:,1)));
[rowS,colS] = find(E_S(:,:,1)==minMatrix);
Swift_Stress = Swift_hardening(Strain,S(rowS,colS,:));
Swift_neck = Swift_necking(S(rowS,colS,:));
SwiftFit(1,:)=S(rowS(1,1),colS(1,1),:);

minMatrix = min(min(E_V(:,:,1)));
[rowV,colV] = find(E_V(:,:,1)==minMatrix);
Voce_Stress = Voce_hardening(Strain,V(rowV,colV,:));
Voce_neck  = Voce_necking(V(rowV,colV,:));
VoceFit(1,:)=V(rowV(1,1),colV(1,1),:);
 
Swift_neck_error = abs(Swift_neck/Strain(end)-1);
Voce_neck_error  = abs(Voce_neck/Strain(end)-1);
    
disp(['Experimental necking strain: ' num2str(Strain(end))]);
disp(['Predicted necking strain (Swift): ' num2str(Swift_neck)]);
disp(['Predicted necking strain (Voce): ' num2str(Voce_neck)]);
disp(['Error (Swift): ' num2str(Swift_neck_error)]);
disp(['Error (Voce): ' num2str(Voce_neck_error)]);
disp([' ']);
SValpha_ini=Voce_neck_error/(Swift_neck_error+Voce_neck_error);
%-------------------------------------------------------
% Plot SV Hardening
%-------------------------------------------------------   
figure(6);

line(Strain,Stress,'LineStyle','none','Color','k','linewidth',1.5,'Marker','.','MarkerSize',10,'MarkerEdgeColor','k'); hold on;
Strain = linspace(Strain(1),1.5*Strain(end),100)';
line(Strain,Swift_hardening(Strain,S(rowS(1,1),colS(1,1),:)),'LineStyle','-','Color','r','linewidth',1.5); hold on;
text(Strain(end),Swift_hardening(Strain(end),S(rowS(1,1),colS(1,1),:)),'\leftarrow Swift'); hold on;

line(Strain,Voce_hardening(Strain,V(rowV(1,1),colV(1,1),:)),'LineStyle','-','Color','g','linewidth',1.5); hold on; 
text(Strain(end),Voce_hardening(Strain(end),V(rowV(1,1),colV(1,1),:)),'\leftarrow Voce'); hold on;

set(gca,'XMinorTick','on');
% set(gca,'XTick',[0:0.05:90]);
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
xmax=1.2*max(Strain);
ymin=0.8*min(Stress);
ymax=1.2*max(Stress);
eval(['axis([0 ' num2str(xmax) ' ' num2str(ymin) ' ' num2str(ymax) '])']);
xlabel('Plastic Strain [-]','FontSize',18);
ylabel('True Stress [MPa]','FontSize',18);
axis square;
box on

clear row* col* S V E_S E_V Swift_* Voce_*
cd _results
savefig(6,'__Swift_Voce_fit');
saveas(6,['__Swift_Voce_fit.png'],'png');
cd ..
%----------------------------------------------------------------------
% MAKE INITIAL GUESS FOR MATERIAL MODEL
%----------------------------------------------------------------------
% %-----------------------------------
% % MAKE INITIAL GUESS FOR HILL48
% %-----------------------------------
% % %-----------------------------------
% % % Calculate P's initial guess
% % %-----------------------------------
% SigY_mean_EB=SigY_mean.UT_00;
% P12_ini = 0.5.*(1.-(SigY_mean.UT_00./SigY_mean_EB).^2-(SigY_mean.UT_00./SigY_mean.UT_90).^2);
% P22_ini = (SigY_mean.UT_00./SigY_mean.UT_90).^2;
% P44_ini = (2.*SigY_mean.UT_00./SigY_mean.UT_45).^2-(SigY_mean.UT_00./SigY_mean_EB).^2;
% %-----------------------------------
% % Calculate G's initial guess
% %-----------------------------------
% G12_ini = -Lank_mean.UT_00./(1+Lank_mean.UT_00);
% G22_ini = Lank_mean.UT_00./Lank_mean.UT_90.*(1+Lank_mean.UT_90)./(1+Lank_mean.UT_00);
% G44_ini = (1.+2.*Lank_mean.UT_45)./Lank_mean.UT_90.*(Lank_mean.UT_00+Lank_mean.UT_90)./(1+Lank_mean.UT_00);
%-----------------------------------
% Calculate P's initial guess only 0-Direction
%-----------------------------------
SigY_mean_EB=SigY_mean.UT_00;
P12_ini = 0.5.*(1.-(SigY_mean.UT_00./SigY_mean_EB).^2-(SigY_mean.UT_00./SigY_mean.UT_00).^2);
P22_ini = (SigY_mean.UT_00./SigY_mean.UT_00).^2;
P44_ini = (2.*SigY_mean.UT_00./SigY_mean.UT_00).^2-(SigY_mean.UT_00./SigY_mean_EB).^2;
%-----------------------------------
% Calculate G's initial guess
%-----------------------------------
G12_ini = -Lank_mean.UT_00./(1+Lank_mean.UT_00);
G22_ini = Lank_mean.UT_00./Lank_mean.UT_00.*(1+Lank_mean.UT_00)./(1+Lank_mean.UT_00);
G44_ini = (1.+2.*Lank_mean.UT_00)./Lank_mean.UT_00.*(Lank_mean.UT_00+Lank_mean.UT_00)./(1+Lank_mean.UT_00);
%-----------------------------------
% Get available Orientations
%-----------------------------------
for i=1:numel(names);
name=char(names(i))   
geometry= name(1:2);
orientation= name(4:5);
adeg(i)=str2num(orientation);
end
%-----------------------------------
% Optimize P's initial guess
%-----------------------------------
global Gi Pi adeg
Pi=cell2mat(struct2cell(SigYRatio));
P_ini=abs([0.5 1 3.]);
options = optimset('Display','off','TolFun',1e-16, 'TolX',1e-16,'MaxFunEvals',1000000,'MaxIter',100000);
% [Gs,fval,exitflag,output] = fminsearchbnd(@CalculateErrorHill, G_ini, lb, ub, options)
[Ps,fval,exitflag,output] = fminsearch(@CalculateErrorHillPs, P_ini, options);

% P12=-Ps(1);
% P22=Ps(2);
% P44=Ps(3);
P12=-0.5;
P22=1;
P44=3.0;
clear P_ini Ps fval exitflag Pi options output
%-----------------------------------
% Optimize G's initial guess
%-----------------------------------
Gi=cell2mat(struct2cell(Lank_mean));
G_ini=abs([0.5 1 3.]);

options = optimset('Display','off','TolFun',1e-16, 'TolX',1e-16,'MaxFunEvals',1000000,'MaxIter',100000);
% [Gs,fval,exitflag,output] = fminsearchbnd(@CalculateErrorHill, G_ini, lb, ub, options)
[Gs,fval,exitflag,output] = fminsearch(@CalculateErrorHillGs, G_ini, options);

% G12=-Gs(1);
% G22=Gs(2);
% G44=Gs(3);
G12=-0.5;
G22=1;
G44=3.;
clear G_ini Gs fval exitflag Gi options output
%-----------------------------------
% Plot Hill 48 SigYratios and r-values
%-----------------------------------
adeg = [0:1:90];
Y_ini=1./(sqrt(cos(adeg.*pi()./180).^4+sin(adeg.*pi()./180).^2.*cos(adeg.*pi()./180).^2.*(2.*P12_ini+P44_ini)+sin(adeg.*pi()./180).^4.*P22_ini));
r_ini=(((G44_ini+2.*G12_ini-G22_ini-1).*sin(adeg.*pi()./180).^2.*cos(adeg.*pi()./180).^2-G12_ini)./((1-G22_ini).*cos(adeg.*pi()./180).^2+G12_ini+G22_ini));
Y=1./(sqrt(cos(adeg.*pi()./180).^4+sin(adeg.*pi()./180).^2.*cos(adeg.*pi()./180).^2.*(2.*P12+P44)+sin(adeg.*pi()./180).^4.*P22));
r=(((G44+2.*G12-G22-1).*sin(adeg.*pi()./180).^2.*cos(adeg.*pi()./180).^2-G12)./((1-G22).*cos(adeg.*pi()./180).^2+G12+G22));
cd _results
fig=openfig('__Y-R_ratios_Angle.fig');

allaxes = findall(fig, 'type', 'axes');
if (strncmp(get(allaxes(1),'YAxisLocation'),'left',4)==1) 
    set(fig,'CurrentAxes',allaxes(1));
    plot(adeg,Y_ini,'--k','linewidth',1.0);
    plot(adeg,Y,'-k','linewidth',1.0);
    set(fig,'CurrentAxes',allaxes(2));
    plot(adeg,r_ini,'--b','linewidth',1.0);
    plot(adeg,r,'-b','linewidth',1.0);
elseif (strncmp(get(allaxes(1),'YAxisLocation'),'right',5)==1)
    set(fig,'CurrentAxes',allaxes(2));
    plot(adeg,Y_ini,'--k','linewidth',1.0);
    plot(adeg,Y,'-k','linewidth',1.0);
    set(fig,'CurrentAxes',allaxes(1));
    plot(adeg,r_ini,'--b','linewidth',1.0);
    plot(adeg,r,'-b','linewidth',1.0);
end
annotation('textbox',...
    [0.2 0.85 0.15 0.05],...
    'String',{'Hill''48'},...
    'FontSize',18,...
    'FontName','Tahoma',...
    'Color',[0 0 0]);

annotation('textbox',...
    [0.2 0.80 0.15 0.05],...
    'EdgeColor','none',...
    'String',{'-- ini'},...
    'FontSize',18,...
    'FontName','Tahoma',...
    'Color',[0 0 0]);
annotation('textbox',...
    [0.2 0.75 0.15 0.05],...
    'EdgeColor','none',...
    'String',{'- opti'},...
    'FontSize',18,...
    'FontName','Tahoma',...
    'Color',[0 0 0]);


savefig(fig,'__Y-R_ratios_Angle_Hill48')
saveas(fig,['__Y-R_ratios_Angle_Hill48.png'],'png');

cd ..

clear Y r adeg 
%-------------------------------------------------------
% SAVE FOR SIMULATION HILL48
%------------------------------------------------------- 
for i=1:numel(names);
name=char(names(i))   
geometry= name(1:2);
orientation= name(4:5);
angle(i)=str2num(orientation);
end
% orientation
output = zeros(4,8);
for i =1:length(angle)
  output(1,i)=angle(i);
end
% LANKFORD
temp=cell2mat(struct2cell(Lank_mean))';
for i=1:numel(temp)
  output(4,i)=temp(i);
end
clear temp
% YIELD STRESSES
temp=cell2mat(struct2cell(SigY_mean))';
for i=1:numel(temp)
  output(2,i)=temp(i);
end
clear temp
% YIELD RATIOS
temp=cell2mat(struct2cell(SigYRatio))';
for i=1:numel(temp)
  output(3,i)=temp(i);
end
clear temp
% INITIAL Ps and Gs
output(5,1:3)=[P12_ini P22_ini P44_ini];
output(5,4:6)=[G12_ini G22_ini G44_ini];
% OPTIMIZED Ps and Gs
output(6,1:3)=[P12 P22 P44];
output(6,4:6)=[G12 G22 G44];
% SWIFT/VOCE FIT
output(7,1:3)=SwiftFit;
output(7,4:6)=VoceFit;
output(7,7)=SValpha_ini;

% disp(output)

mkdir _simulation
cd _simulation
dlmwrite(['_Hill48naFr_pre.csv'],[output],'precision','%.6f');
cd ..

cd _results

dlmwrite(['_Hill48naFr_pre.csv'],[output],'precision','%.6f');
cd ..
%-----------------------------------
% MAKE INITIAL GUESS FOR YLD2000
%-----------------------------------
global a0 P Pi Sigma K0;

Pi=[SigYRatio.UT_00; Lank_mean.UT_00; ...
    SigYRatio.UT_00; Lank_mean.UT_00; ...
   SigYRatio.UT_00; Lank_mean.UT_00; ...
    SigYRatio.UT_00; Lank_mean.UT_00; ...
    SigYRatio.UT_00; Lank_mean.UT_00; ...
    SigYRatio.UT_00; Lank_mean.UT_00; ...
    SigYRatio.UT_00; Lank_mean.UT_00; ...
    YEBT; rEBT]; %  

%a0 = a0. %exponent

[alpha_pre,P_f,error_f] = AluminumCalibration(Pi,a0)
% alpha_pre=[0.940631008444053,0.998708806643439,0.890552926298443,1.0059885616344,0.993486594914199, 0.880368970241785, 0.955198047310224, 1.1671987325567];
%-----------------------------------
% PLOT 3D YLD2000
%-----------------------------------
phi   = linspace(0,360,180);
theta = linspace(0,90,45);

[Phi, Theta] = meshgrid(phi/180*pi,theta/180*pi);
X = sin(Theta).*cos(Phi);
Y = sin(Theta).*sin(Phi);
Z = cos(Theta).*ones(size(Phi));

Sigma_all= Yld2000_plot3d(X,Y,Z,alpha_pre,a0);
Scale = Sigma_all(end,end);

figure(7);
axeb=axes('FontSize',18);

surf(Sigma_all.*X/Scale,Sigma_all.*Y/Scale,Sigma_all.*Z/Scale, 'LineStyle','none','FaceAlpha', 0.8); hold on

X2D = Sigma_all(end,:).*X(end,:)/Scale;
Y2D = Sigma_all(end,:).*Y(end,:)/Scale;

view([116,33]);
set(axeb,'XLim',[-1.5 1.5]);
set(axeb,'YLim',[-1.5 1.5]);
set(axeb,'ZLim',[0 0.75]);
set(axeb,'XTick',[-1.5:0.5:1.5]);
set(axeb,'XTicklabel', num2str(get(axeb, 'XTick')', '%.1f'));
set(axeb,'YTick',[-1.5:0.5:1.5]);
set(axeb,'YTicklabel', num2str(get(axeb, 'YTick')', '%.1f'));
set(axeb,'ZTick',[0.:0.25:0.75]);
set(axeb,'ZTicklabel', num2str(get(axeb, 'ZTick')', '%.2f'));

xlabel(axeb,'\sigma_x_x/\sigma_y_,_0 [-]'); 
ylabel(axeb,'\sigma_y_y/\sigma_y_,_0 [-]'); 
zlabel(axeb,'\sigma_x_y/\sigma_y_,_0 [-]'); 
box(axeb,'on');
hold(axeb,'all');

%-----------------------------------
% PLOT 2D YLD2000 INI
%-----------------------------------
figure(5);
axeb=axes('FontSize',18);

line(X2D,Y2D,'Color','k'); hold on

set(axeb,'XLim',[-1.5 1.5]);
set(axeb,'YLim',[-1.5 1.5]);
set(axeb,'XTick',[-1.5:0.5:1.5]);
set(axeb,'XTicklabel', num2str(get(axeb, 'XTick')', '%.1f'));
set(axeb,'YTick',[-1.5:0.5:1.5]);
set(axeb,'YTicklabel', num2str(get(axeb, 'YTick')', '%.1f'));
xlabel(axeb,'\sigma_x_x/\sigma_y_,_0 [-]'); 
ylabel(axeb,'\sigma_y_y/\sigma_y_,_0 [-]');  
axis square
grid on
box(axeb,'on');
hold(axeb,'all');

%-----------------------------------
% PLOT 2D YLD2000 INI in SigY/R plot
%-----------------------------------
[Mat_Ori,Y,r]=Yld20002d_yld_r_change(alpha_pre,a0);
figure(2);
set(2,'CurrentAxes',ax1);
plot(Mat_Ori*180/pi,Y,'k','linewidth',1.0);
set(2,'CurrentAxes',ax2);
plot(Mat_Ori*180/pi,r,'b','linewidth',1.0);

annotation('textbox',...
    [0.2 0.85 0.20 0.05],...
    'String',{'YLD2000-3D'},...
    'FontSize',18,...
    'FontName','Tahoma',...
    'Color',[0 0 0]);

%-------------------------------------------------------
% SAVE FOR SIMULATION 
%------------------------------------------------------- 
for i=1:numel(names);
name=char(names(i))   
geometry= name(1:2);
orientation= name(4:5);
angle(i)=str2num(orientation);
end
% orientation
output = zeros(4,8);
for i =1:length(angle)
  output(1,i)=angle(i);
end
% LANKFORD
temp=cell2mat(struct2cell(Lank_mean))';
for i=1:numel(temp)
  output(4,i)=temp(i);
end
clear temp
% YIELD STRESSES
temp=cell2mat(struct2cell(SigY_mean))';
for i=1:numel(temp)
  output(2,i)=temp(i);
end
clear temp
% YIELD RATIOS
temp=cell2mat(struct2cell(SigYRatio))';
for i=1:numel(temp)
  output(3,i)=temp(i);
end
clear temp
% INITIAL Ps and Gs
for i =1:length(alpha_pre)
  output(5,i)=alpha_pre(1,i);
end
% OPTIMIZED Ps and Gs
output(6,1:6)=[0 0 0 0 0 0];
% SWIFT/VOCE FIT
output(7,1:3)=SwiftFit;
output(7,4:6)=VoceFit;
output(7,7)=SValpha_ini;

disp(output)

mkdir _simulation
cd _simulation
dlmwrite(['_YLD2000_pre.csv'],[output],'precision','%.6f');
cd ..

cd _results

dlmwrite(['_YLD2000_pre.csv'],[output],'precision','%.6f');

savefig(2,'__Y-R_ratios_Angle_YLD2000')
saveas(2,['__Y-R_ratios_Angle_YLD2000.png'],'png');
savefig(7,'__YLD2000_plot_3D');
saveas(7,['__YLD2000_plot_3D.png'],'png');
savefig(5,'__YLD2000_3D');
saveas(5,['__YLD2000_3D.png'],'png');

cd ..
 end

%-------------------------------------------------------
% Auxiliary Functions
%------------------------------------------------------- 

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

function error=CalculateErrorHillGs(G_ini)

    global Gi Pi adeg
 
    G12 = -G_ini(1);
    G22 = G_ini(2);
    G44 = G_ini(3);
    
    a=adeg.*pi()./180;
    
    r=((G44+2.*G12-G22-1.).*sin(a).^2.*cos(a).^2-G12)./ ...
        ((1-G22).*cos(a).^2+G12+G22);
      
    for i=1:numel(adeg)
%     for i=[2 3 4 5 6]
    error(i)=abs((r(i)./Gi(i))-1).^2;%.^2;
    end
    
    error=sum(error);
       
end

function error=CalculateErrorHillPs(P_ini)

    global Gi Pi adeg
 
    P12 = -P_ini(1);
    P22 = P_ini(2);
    P44 = P_ini(3);
    
    a=adeg.*pi()./180;
    
    sig_r=1./(sqrt(cos(adeg.*pi()./180).^4+sin(adeg.*pi()./180).^2.* ...
        cos(adeg.*pi()./180).^2.*(2.*P12+P44)+sin(adeg.*pi()./180).^4.*P22));
      
    for i=1:numel(adeg)
%     for i=[2 3 4 5 6]
    error(i)=abs((sig_r(i)./Pi(i))-1).^2;%.^2;
    end
    
    error=sum(error);
       
end
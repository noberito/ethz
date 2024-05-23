function Hill48_preprocessing
clear all; close all; clc; format long
minPStrain=0.002;
%-----------------------------------
% Load result filenames
%-----------------------------------
cd _results
listOf = dir('*.csv');
for i=1:numel(listOf)
filenames{i}=listOf(i).name;
end
%-----------------------------------
count=1; 
for i=1:numel(filenames)   
 
if (strncmp(char(filenames(i)),'UT00',4)==1)
    
    A=csvread(char(filenames(i)),1,0);
    Lankford(count).deg00=A(1,13);   
    
    a=A(:,10);
    aa=any(a~=0,2);
    
    PlasticStrain(count).deg00=A(aa,9);
    PlasticStress(count).deg00=A(aa,10);
    PlasticStrain_max(i)=max(PlasticStrain(count).deg00);
    count=count+1; 
    clear A a aa
end
end
count=1; 
for i=1:numel(filenames)  
if (strncmp(char(filenames(i)),'UT45',4)==1)  
    
    A=csvread(char(filenames(i)),1,0);
    Lankford(count).deg45=A(1,13);   
    
    a=A(:,10);
    aa=any(a~=0,2);
    
    PlasticStrain(count).deg45=A(aa,9);
    PlasticStress(count).deg45=A(aa,10);
    PlasticStrain_max(i)=max(PlasticStrain(count).deg45);
    count=count+1; 
    clear A a aa    
end
end
count=1; 
for i=1:numel(filenames)  
if (strncmp(char(filenames(i)),'UT90',4)==1)  
    
    A=csvread(char(filenames(i)),1,0);
    Lankford(count).deg90=A(1,13);   
    
    a=A(:,10);
    aa=any(a~=0,2);
    
    PlasticStrain(count).deg90=A(aa,9);
    PlasticStress(count).deg90=A(aa,10);
    PlasticStrain_max(i)=max(PlasticStrain(count).deg90);
    count=count+1; 
    clear A a aa
    
end
end
clear count 
cd ..
%-----------------------------------
% Find smallest of mac plastic strains for spline 
% TO BE IMPROVED WITH STATISTICAL METHODS
%-----------------------------------
% if min(PlasticStrain_max)>0.1
%     PlasticStrain_int=linspace(minPStrain,0.1,100)';    
% else
    PlasticStrain_int=(minPStrain:0.001:.1)';
%     PlasticStrain_int=linspace(minPStrain,min(PlasticStrain_max),100)';
% end

[zz j] = size(PlasticStrain);
for i=1:j
    PlasticStress_int(i).deg00=spline(PlasticStrain(i).deg00,PlasticStress(i).deg00,PlasticStrain_int);
    PlasticStress_int(i).deg45=spline(PlasticStrain(i).deg45,PlasticStress(i).deg45,PlasticStrain_int);
    PlasticStress_int(i).deg90=spline(PlasticStrain(i).deg90,PlasticStress(i).deg90,PlasticStrain_int);
    SigY(i).deg00 = PlasticStress_int(i).deg00;
    SigY(i).deg45 = PlasticStress_int(i).deg45;
    SigY(i).deg90 = PlasticStress_int(i).deg90;
end
SigY_00=0;
SigY_45=0;
SigY_90=0;
r_00=0;
r_45=0;
r_90=0;
% WHERE TO CALC SIG RATIOS
ratiolocation =2 
for i=1:j
plot(PlasticStrain_int,SigY(i).deg00./SigY(i).deg00);hold on
pos=randi([10,length(PlasticStrain_int)-10],1);
eval(['text(PlasticStrain_int(pos),SigY(i).deg00(pos)./SigY(i).deg00(pos)+0.002,''\downarrow UT00_' num2str(i) ' '')']);
plot(PlasticStrain_int,SigY(i).deg45./SigY(i).deg00);
pos=randi([10,length(PlasticStrain_int)-10],1);
eval(['text(PlasticStrain_int(pos),SigY(i).deg45(pos)./SigY(i).deg00(pos)+0.002,''\downarrow UT45_' num2str(i) ' '')']);
plot(PlasticStrain_int,SigY(i).deg90./SigY(i).deg00);
pos=randi([10,length(PlasticStrain_int)-10],1);
eval(['text(PlasticStrain_int(pos),SigY(i).deg90(pos)./SigY(i).deg00(pos)+0.002,''\downarrow UT90_' num2str(i) ' '')']);
plot(PlasticStrain_int(ratiolocation),[0.9:0.01:1.1],'.r');
%SET IMAGE PROPERTIES
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
axis([0 0.1 0.95 1.2]);
xlabel('Plastic Strain [-]','FontSize',18);
ylabel('\sigma_0_,_4_5_,_9_0/\sigma_0 [-]','FontSize',18);
axis square;
box on

SigY_00=SigY_00+SigY(i).deg00(ratiolocation);
SigY_45=SigY_45+SigY(i).deg45(ratiolocation);
SigY_90=SigY_90+SigY(i).deg90(ratiolocation);
r_00=r_00+Lankford(i).deg00;
r_45=r_45+Lankford(i).deg45;
r_90=r_90+Lankford(i).deg90;
end

SigY_00=SigY_00/j
% SigY_45=SigY_45/j
SigY_45=SigY(1).deg45(ratiolocation)
SigY_90=SigY_90/j
r_00=r_00./j
% r_45=r_45./j;
r_45=Lankford(1).deg45
r_90=r_90./j

% TEST
% plot(PlasticStrain_int,SigY(i).deg00);hold on
% plot(PlasticStrain_int,SigY(i).deg45);
% plot(PlasticStrain_int,SigY(i).deg90);
%-----------------------------------
% Calculate P's and G's
%-----------------------------------
P12 = 0.5.*(1.-(SigY_00./SigY_90).^2-1)
P22 = (SigY_00./SigY_90).^2
P44 = (2.*SigY_00./SigY_45).^2-1
% 
G12 = -r_00./(1+r_00)
G22 = r_00./r_90.*(1+r_90)./(1+r_00)
G44 = (1.+2.*r_45)./r_90.*(r_00+r_90)./(1+r_00)

%-------------------------------------------------------
% Start fitting hardening coefficients
%-------------------------------------------------------    
%     
Strain = PlasticStrain(1).deg00;
Stress = PlasticStress(1).deg00;
    
[S,V,E_S,E_V]=FitHardening(Strain,Stress);

minMatrix = min(min(E_S(:,:,1)));
[rowS,colS] = find(E_S(:,:,1)==minMatrix);
Swift_Stress = Swift_hardening(Strain,S(rowS,colS,:));
Swift_neck = Swift_necking(S(rowS,colS,:));
SwiftFit(1,:)=S(rowS,colS,:);

minMatrix = min(min(E_V(:,:,1)));
[rowV,colV] = find(E_V(:,:,1)==minMatrix);
Voce_Stress = Voce_hardening(Strain,V(rowV,colV,:));
Voce_neck  = Voce_necking(V(rowV,colV,:));
VoceFit(1,:)=V(rowV,colV,:);
 
Swift_neck_error = abs(Swift_neck/Strain(end)-1);
Voce_neck_error  = abs(Voce_neck/Strain(end)-1);
    
disp(['Experimental necking strain: ' num2str(Strain(end))]);
disp(['Predicted necking strain (Swift): ' num2str(Swift_neck)]);
disp(['Predicted necking strain (Voce): ' num2str(Voce_neck)]);
disp(['Error (Swift): ' num2str(Swift_neck_error)]);
disp(['Error (Voce): ' num2str(Voce_neck_error)]);
disp([' ']);
alpha_ini=Voce_neck_error/(Swift_neck_error+Voce_neck_error);

%% Plot    
figure1 = figure('Name','Swift/Voce Fit');
axes1 = axes('Parent',figure1,...
    'XMinorTick','on',...
    'YMinorTick','on',...
    'TickLength',[0.01 0.01],...
    'LineWidth',1,...
    'FontSize',18,...
    'XColor','k',...
    'YColor','k');
%     'XTick',[0:0.2:1.1],...
%     'YTick',[0:15:90],...
xlabel('Plastic Strain [-]','FontSize',18);
ylabel('True Stress [MPa]','FontSize',18);

axis square;
box on;

line(Strain,Stress,'LineStyle','none','Color','k','linewidth',1.5,'Marker','.','MarkerSize',10,'MarkerEdgeColor','k'); hold on;
Strain = linspace(Strain(1),1.5*Strain(end),100)';
line(Strain,Swift_hardening(Strain,S(rowS,colS,:)),'LineStyle','-','Color','r','linewidth',1.5); hold on;
text(Strain(end),Swift_hardening(Strain(end),S(rowS,colS,:)),'\leftarrow Swift'); hold on;

line(Strain,Voce_hardening(Strain,V(rowV,colV,:)),'LineStyle','-','Color','g','linewidth',1.5); hold on; 
text(Strain(end),Voce_hardening(Strain(end),V(rowV,colV,:)),'\leftarrow Voce'); hold on;
% eval(['axis([0 ' num2str(1.1*Strain(end)) num2str(0.9*Stress(1)) num2str(1.1*Stress(end)) '])']);

fileName=['SV_Fit.tif'];
print('-dtiff','-r220',fileName);

mkdir _simulation
cd _simulation

dlmwrite(['Hill48naFr_SV_Ini.csv'],[SigY_00,SigY_45,SigY_90,r_00,r_45,r_90,0;...
    SwiftFit,VoceFit,alpha_ini;...
    P12,P22,P44,G12,G22,G44,0],'precision','%.6f');
cd ..
cd _results
dlmwrite(['Hill48naFr_SV_Ini.csv'],[SigY_00,SigY_45,SigY_90,r_00,r_45,r_90,0;...
    SwiftFit,VoceFit,alpha_ini;...
    P12,P22,P44,G12,G22,G44,0],'precision','%.6f');
savefig(1,'10_Y_ratios')
saveas(1,['10_Y_ratios.png'],'png');
savefig(2,'11_Swift_Voce_fit');
saveas(2,['11_Swift_Voce_fit.png'],'png');

cd ..
% disp(S);
% disp(V);
 end
% 
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

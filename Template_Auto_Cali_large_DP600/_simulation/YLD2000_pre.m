function YLD2000_pre
global a0 P Pi Sigma K0 alpha_pre;
clear all; close all; clc; format long
minPStrain=0.002;
PercentStrain=0.03;
YieldStrain = 0.002;
%-----------------------------------
% Load results from files
%-----------------------------------
cd _results
listOf = dir('*.csv');
ii=1;
for i=1:numel(listOf)
if (strncmp(listOf(i).name,'UT',2)==1)
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
    
    eval(['PlasticStrain_int_max_all = min(PlasticStrain_int_max_all,max(PlasticStrain_int(testnumber).' nameori '));']); 
        
    clear A a aa pos name geometry orientation plotname testname testnumber nameori

end
cd ..
%
PlasticStrain_int_all=(minPStrain:0.0001:PlasticStrain_int_max_all)'
ratiolocation =find(PlasticStrain_int_all==PercentStrain);
yieldlocation =find(PlasticStrain_int_all==YieldStrain);
%
names=fieldnames(PlasticStress_int);
% SPLINE ALL SIG?EPSILON UP TO MAX EPS OF ALL
for i=1:numel(names);
name=char(names(i))   
geometry= name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];
for testnumber=1:100
try
eval(['PlasticStress_int_all(testnumber).' nameori '=pchip(PlasticStrain(testnumber).' nameori ',PlasticStress(testnumber).' nameori ',PlasticStrain_int_all);']);
end
end
clear i A a aa pos name geometry orientation plotname testname testnumber nameori
end
%-----------------------------------
% MEAN OF STRESS & YIELD STRESS
%-----------------------------------
for i=1:numel(names);
name=char(names(i))   
geometry=  name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];
for ii=1:100
try
eval(['temp=PlasticStress_int(ii).' nameori ';']);  
eval(['SigY(ii).' nameori '=PlasticStress_int(ii).' nameori '(ratiolocation);']);
if temp>0
ll=ii;
end
end
clear temp

eval(['PlasticStress_int_mean.' nameori '=mean(cat(ll,PlasticStress_int_all.' nameori '),2);']);

end
clear ll i ii name geometry orientation nameori
end
%-----------------------------------
% PLOT SIG RATIOS
%-----------------------------------
figure(1);
for i=1:numel(names);
name=char(names(i))   
geometry=  name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];
plotname =[geometry '-' orientation];
eval(['SigYRatio_mean.' nameori '=PlasticStress_int_mean.' nameori './PlasticStress_int_mean.UT_00;']);
eval(['SigYRatio.' nameori '=SigYRatio_mean.' nameori '(ratiolocation);']);
eval(['plot(PlasticStrain_int_all,SigYRatio_mean.' nameori ');']); hold on
pos=randi([10,length(PlasticStrain_int_all)-10],1);
eval(['text(PlasticStrain_int_all(pos),mean(PlasticStress_int_mean.' nameori '(pos-10:pos+10)./PlasticStress_int_mean.UT_00(pos-10:pos+10)),'' ' plotname ' '')']);
end
plot(PlasticStrain_int_all(ratiolocation),[0.9:0.01:1.1],'.r');
%SET IMAGE PROPERTIES
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
xmax=max(PlasticStrain_int_all);
ymin=min(min(cell2mat(struct2cell(SigYRatio_mean))))-0.05;
ymax=max(max(cell2mat(struct2cell(SigYRatio_mean))))+0.05;
eval(['axis([0 ' num2str(xmax) '  ' num2str(ymin) ' ' num2str(ymax) '])']);
xlabel('Plastic Strain [-]','FontSize',18);
ylabel('\sigma_0_,_4_5_,_9_0/\sigma_0 [-]','FontSize',18);
axis square;
box on
%-----------------------------------
% MEAN OF LANKFORD AND SIG
%-----------------------------------
names=fieldnames(PlasticStress_int);
for i=1:numel(names);
name=char(names(i))   
geometry=  name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];

eval(['SigY_mean.' nameori '=mean([SigY.' nameori ']);']);   
eval(['SigYRatio.' nameori '=mean([SigYRatio.' nameori ']);']);  
eval(['Lank_mean.' nameori '=mean([Lankford.' nameori ']);']);

for ii=1:100
try
eval(['temp=PlasticStress_int_all(ii).' nameori ';']);  
if temp>0
ll=ii;
end
end
eval(['PlasticStress_int_mean.' nameori '=mean(cat(ll,PlasticStress_int_all.' nameori '),2);']);

end
clear ll ii name geometry orientation nameori temp
end
figure(2);
for i=1:numel(fieldnames(PlasticStress))
name=char(names(i))   
geometry=  name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];
plotname =[geometry '-' orientation];

    eval(['plot(str2num(orientation),SigYRatio.' nameori ',''.k'');']); hold on
%PLOT MAX/MIN
    eval(['plot(str2num(orientation),max(SigYRatio_mean.' nameori '),''xk'');']); hold on
    eval(['plot(str2num(orientation),min(SigYRatio_mean.' nameori '),''xk'');']); hold on
%PLOT ALL
%     eval(['plot(str2num(orientation),(SigYRatio_mean.' nameori '),''xk'');']); hold on
    clear A a aa pos
    
end
set(gca,'XMinorTick','on');
set(gca,'XTick',[0:15:90]);
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
ymin=0.9*min(cell2mat(struct2cell(SigYRatio)));
ymax=1.1*max(cell2mat(struct2cell(SigYRatio)));
eval(['axis([-5 95 ' num2str(ymin) ' ' num2str(ymax) '])']);
xlabel('Material Orientation [°]','FontSize',18);
ylabel('\sigma_x_°/\sigma_0_° [-]','FontSize',18);
axis square;
box on
%-----------------------------------
% Plot R-Values
%-----------------------------------
figure(3);
for i=1:numel(fieldnames(Lankford))
name=char(names(i))   
geometry=  name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];
plotname =[geometry '-' orientation];    

    eval(['plot(str2num(orientation),Lank_mean.' nameori ',''.k'');']); hold on
    
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
eval(['plot(str2num(orientation),(cat(ll,Lankford.' nameori ')),''xk'');']); hold on

clear ll ii name geometry orientation nameori temp
    
end
set(gca,'XMinorTick','on');
set(gca,'XTick',[0:15:90]);
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
ymin=0.9*min(min(cell2mat(struct2cell(Lankford))));
ymax=1.1*max(max(cell2mat(struct2cell(Lankford))));
eval(['axis([-5 95 ' num2str(ymin) ' ' num2str(ymax) '])']);
xlabel('Material Orientation [°]','FontSize',18);
ylabel('r-Value [-]','FontSize',18);
axis square;
box on
clear plotname i xmax ymax ymin 
% %-----------------------------------
% % CALCULATE AVERAGE SIGYRATIOS AND R-VALUES????
% %-----------------------------------
for i=1:numel(fieldnames(SigY))
nameori=char(names(i));
angle(i)=str2num(nameori(4:5));


 eval(['SigY_pre(i)=SigY_mean.' nameori ';']); hold on
 eval(['Lank_pre(i)=Lank_mean.' nameori ';']); hold on
 
end
clear i
Lank_pre(1)=Lankford(1).UT_00; %CORRECTION FOR BAD UT00_2
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

    eval(['A(:,1)=PlasticStrain(1).' nameori ';']);
    eval(['A(:,2)=PlasticStress(1).' nameori ';']);

name = strcat(nameori,'_exp.csv')
dlmwrite(name,A);
clear A
end
cd ..
%-----------------------------------
% MAKE INITIAL GUESS FOR YLD2000
%-----------------------------------

% Pi=[SigY(1).UT_00; Lankford(1).UT_00; ...
%     SigY(1).UT_15; Lankford(1).UT_15; ...
%     SigY(1).UT_30; Lankford(1).UT_30; ...
%     SigY(1).UT_45; Lankford(1).UT_45; ...
%     SigY(1).UT_60; Lankford(1).UT_60; ...
%     SigY(1).UT_75; Lankford(1).UT_75; ...
%     SigY(1).UT_90; Lankford(1).UT_90]; %SigY(1).UT_00; 1.; ...

Pi=[SigY_mean.UT_00; Lankford(1).UT_00; ...
    SigY_mean.UT_45; Lank_mean.UT_45; ...
    SigY_mean.UT_90; Lankford(1).UT_90; ...
    SigY_mean.UT_00; 1.]; 

a0 = 8 %exponent

[alpha_pre,P_f,error_f] = AluminumCalibration(Pi,a0)

% P=zeros(length(Pi),1);
% K0=Pi(1);
% a=[1 1 1 1 1 1 1 1 1 1];
 
% Sigma=[1/2;1/2;1/2];
% [Phi,dPhi]=EvalPhi(100,a);
% options = optimset('Display','off','TolFun',1e-12, 'TolX',1e-12,'MaxFunEvals',10000,'MaxIter',10000);
% [A,fval,exitflag,output] = fminsearch(@CalculateError,a,options);
% 
% alpha=A;
% P_f=[P Pi P./Pi];
% error=fval;
% error_f=error
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

figure(4);
axeb=axes('FontSize',18);

surf(Sigma_all.*X/Scale,Sigma_all.*Y/Scale,Sigma_all.*Z/Scale, 'LineStyle','none','FaceAlpha', 0.8)

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

line(X2D,Y2D,'Color','k'), hold on

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
figure(3);
plot(Mat_Ori*180/pi,r,'b','linewidth',1.0);
figure(2);
plot(Mat_Ori*180/pi,Y,'k','linewidth',1.0);

%-------------------------------------------------------
% Start fitting hardening coefficients
%-------------------------------------------------------    
% CHOOSE the right UT000 which criterion???
Strain = PlasticStrain(1).UT_00;
Stress = PlasticStress(1).UT_00;
    
Strain = Strain(1:350,:);
Stress = Stress(1:350,:);

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
figure(6);

line(Strain,Stress,'LineStyle','none','Color','k','linewidth',1.5,'Marker','.','MarkerSize',10,'MarkerEdgeColor','k'); hold on;
Strain = linspace(Strain(1),1.5*Strain(end),100)';
line(Strain,Swift_hardening(Strain,S(rowS,colS,:)),'LineStyle','-','Color','r','linewidth',1.5); hold on;
text(Strain(end),Swift_hardening(Strain(end),S(rowS,colS,:)),'\leftarrow Swift'); hold on;

line(Strain,Voce_hardening(Strain,V(rowV,colV,:)),'LineStyle','-','Color','g','linewidth',1.5); hold on; 
text(Strain(end),Voce_hardening(Strain(end),V(rowV,colV,:)),'\leftarrow Voce'); hold on;

set(gca,'XMinorTick','on');
% set(gca,'XTick',[0:0.05:90]);
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
xmax=1.2*max(Strain);
ymin=0.8*min(Stress);
ymax=1.3*max(Stress);
eval(['axis([0 ' num2str(xmax) ' ' num2str(ymin) ' ' num2str(ymax) '])']);
xlabel('Plastic Strain [-]','FontSize',18);
ylabel('True Stress [MPa]','FontSize',18);
axis square;
box on

% eval(['axis([0 ' num2str(1.1*Strain(end)) num2str(0.9*Stress(1)) num2str(1.1*Stress(end)) '])']);

%-------------------------------------------------------
% SAVE DATA AND PLOTS
%-------------------------------------------------------    
output = zeros(4,8);
for i =1:length(angle)
  output(1,i)=angle(1,i);
end
for i =1:length(Lank_pre)
  output(2,i)=Lank_pre(1,i);
end
for i =1:length(alpha_pre)
  output(3,i)=alpha_pre(1,i);
end
output(4,1:3)=SwiftFit;
output(4,4:6)=VoceFit;
output(4,7)=alpha_ini;

mkdir _simulation
cd _simulation
dlmwrite(['YLD2000_pre.csv'],[output],'precision','%.6f');
cd ..

cd _results

dlmwrite(['YLD2000_pre.csv'],[output],'precision','%.6f');

savefig(1,'__Y_ratios_Strain')
saveas(1,['__Y_ratios_Strain.png'],'png');
savefig(2,'__Y_ratios_Angle')
saveas(2,['__Y_ratios_Angle.png'],'png');
savefig(3,'__R-Values_Angle')
saveas(3,['__R-Values_Angle.png'],'png');
savefig(4,'__YLD2000_plot_3D');
saveas(4,['__YLD2000_plot_3D.png'],'png');
savefig(5,'__YLD2000_2D');
saveas(5,['__YLD2000_2D.png'],'png');
savefig(6,'__Swift_Voce_fit');
saveas(6,['__Swift_Voce_fit.png'],'png');

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
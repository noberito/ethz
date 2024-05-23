%function UT_post(filename,Young,Poisson,minPStrain)
clear all; close all; clc;
cd _data
%-----------------------------------
% FIND AND LOAD UT FILES
%-----------------------------------
currentFolder =pwd;
pathparts = strsplit(currentFolder,filesep);
listOf = dir('*.txt');
ii=1;
for i=1:numel(listOf)
str=listOf(i).name;
if (strfind(str,'UT_'))
filenames{ii}=listOf(i).name;
ii=ii+1;
end
end
NbTest = numel(filenames);
clear i ii listOf str

YoungId     ='Theo'; %'Exp' or 'Theo' for use of experimental or theoretical Young modulus
YoungTheo   = 205000; % [MPa]
PoissonTheo = 0.3;
minPStrain  = 0.0005; % Minimum plastic strain use for the hardening calibration

xaxis=zeros(1,20);
yaxis=zeros(1,20);
cd ..
%-----------------------------------
% Identify specimen type 
%-----------------------------------
for zz=1:NbTest
cd _data    
filename=filenames{zz};
%-----------------------------------
% Load data from file
%-----------------------------------
fid=fopen(filename,'r');
AA = textscan(fid,'%s',17,'Delimiter',';');
A  = textscan(fid,'%s',17,'Delimiter',';');
Material=cell2mat(A{1}(2));
geometry=cell2mat(A{1}(3));
orientation=cell2mat(A{1}(4));
testnumber=cell2mat(A{1}(5));
NominalSpeed=str2num(cell2mat(A{1}(6)));
thickness =str2num(cell2mat(A{1}(7)));%2.55;
outerwidth =str2num(cell2mat(A{1}(8))); %10;
innerwidth = str2num(cell2mat(A{1}(9))); 
mmperpix=str2num(cell2mat(A{1}(10)));
%TO BE ADDED FROM MEASUREMENT
extenso_zone_1_4_u1=str2num(cell2mat(A{1}(11)));
extenso_zone_2_3_u2=str2num(cell2mat(A{1}(12)));
extenso_zone_5_6_u3=str2num(cell2mat(A{1}(13)));
fclose(fid);
clear A*;
%
delimiterIn = ';';
headerlinesIn = 4;
A = importdata(filename,delimiterIn,headerlinesIn);
A = A.data;   
cd ..
%-----------------------------------
% TESTNAME & FIND LOSS OF CORRELATION
%-----------------------------------
plotname = [geometry '-' orientation '-' testnumber];
testname = [geometry '_' orientation '_' testnumber];
testname = erase(filename,"_Front_DIC_File.txt");
testname = [cell2mat(pathparts(6)) '_' cell2mat(pathparts(7)) '_' testname '_full'];
for i=1:3
a=A(:,i+1);
aa = isfinite(a);
B =a(aa,:);
bb=any(B~=-1,2);
C=B(bb,:);
lossofcorr(i)=length(C);
clear a aa B bb C
end   
%-----------------------------------
% Force displacement 
%-----------------------------------
if lossofcorr(1)<lossofcorr(2) || lossofcorr(1)<lossofcorr(3)
    lossofcorr(1)=min(lossofcorr(2),lossofcorr(3));
end
Time = A(1:lossofcorr(1),1);
Force = A(1:lossofcorr(1),5);
% Force = smooth(A(1:lossofcorr(1),5),'lowess');
u1=A(1:lossofcorr(1),2);
u2=A(1:lossofcorr(2),3);
u3=A(1:lossofcorr(3),4);
exx= A(1:lossofcorr(1),19);
eyy= A(1:lossofcorr(1),20);
exy= A(1:lossofcorr(1),21);
e1= A(1:lossofcorr(1),22);
e2= A(1:lossofcorr(1),23);
clear A;
figure(1);%=figure('Name','Force Displacement');
plot(u1,Force(1:length(u1)),'.','MarkerSize',4); hold on; 
pos=randi([length(u1)-10,length(u1)-5],1);
eval(['text(u1(pos),Force(pos),''\leftarrow ' plotname ' '')']);
if zz==zz(end)
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
% set(gca,'FontSize',18);
% set(gca,'FontSize',18);
xaxis(1)=max(xaxis(1),1.1*max(u1(:)));
yaxis(1)=max(yaxis(1),1.1*max(Force));
eval(['axis([0 ' num2str(xaxis(1)) ' 0 ' num2str(yaxis(1)) '])']);
% eval(['ylim(axes1,[0 ' num2str(1.1*max(Force)) '])']);
xlabel('Displacement [mm]','FontSize',18);
ylabel('Force [kN]','FontSize',18);
axis square;
box on
end
%PLOT TIME DISPLACEMENT 
figure(7);
plot(Time,u1,'.','MarkerSize',4); hold on
pos=randi([length(Time)-5,length(Time)],1);
eval(['text(Time(pos),u1(pos),''\leftarrow ' plotname ' '')']);
if zz==zz(end)
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
% set(gca,'FontSize',18);
% set(gca,'FontSize',18);
xaxis(7)=max(xaxis(7),1.1*max(Time(end-20:end)));
yaxis(7)=max(yaxis(7),1.1*max(u1(:)));
eval(['axis([0 ' num2str(xaxis(7)) ' 0 ' num2str(yaxis(7)) '])']);
xlabel('Time [s]','FontSize',18);
ylabel('Displacement [mm]','FontSize',18);
axis square;
box on
end
%-----------------------------------
% Eng Stress Strain 
%-----------------------------------
EngStrain_longi=exp(e1)-1;
EngStrain_width=exp(e2)-1;
EngStress=Force(1:length(e1)).*1000./(thickness.*outerwidth);

figure(2)%=figure('Name','Eng. Stress/ Eng. Strain');
plot(EngStrain_longi,EngStress,'.','MarkerSize',4); hold on
pos=randi([length(EngStrain_longi)-5,length(EngStrain_longi)],1);
eval(['text(EngStrain_longi(pos),EngStress(pos),''\leftarrow ' plotname ' '')']);
if zz==zz(end)
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
% set(gca,'FontSize',18);
% set(gca,'FontSize',18);
xaxis(2)=max(xaxis(2),1.1*max(abs(EngStrain_longi((end-2)-20:(end-2)))));
yaxis(2)=max(yaxis(2),1.1*max(EngStress));
eval(['axis([0 ' num2str(xaxis(2)) ' 0 ' num2str(yaxis(2)) '])']);
xlabel('Eng. Strain [-]','FontSize',18);
ylabel('Eng. Stress [MPa]','FontSize',18);
axis square;
box on
end
%-----------------------------------
% True Stress Strain 
%-----------------------------------
[fmax ind]=max(Force);

EngStrain_longi_max=EngStrain_longi(1:min(lossofcorr(2),ind));
EngStrain_width_max=EngStrain_width(1:min(lossofcorr(3),ind));
EngStress_max=EngStress(1:min(lossofcorr(1),ind));
TrueTime=Time(1:min(lossofcorr(1),ind));

TrueStrain_longi=log(1+EngStrain_longi_max);
TrueStrain_width=log(1+EngStrain_width_max);
TrueStress=(1+EngStrain_longi_max).*EngStress_max;

figure(3)%figure3=figure('Name','True Stress/ True Strain');
plot(TrueStrain_longi,TrueStress); hold on
pos=length(TrueStrain_longi);%randi([length(TrueStrain_longi)-50,length(TrueStrain_longi)],1);
eval(['text(TrueStrain_longi(pos),TrueStress(pos),''\leftarrow ' plotname ' '')']);
if zz==zz(end)
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
% set(gca,'FontSize',18);
% set(gca,'FontSize',18);
xaxis(3)=max(xaxis(3),1.1*TrueStrain_longi(end));
yaxis(3)=max(yaxis(3),1.1*TrueStress(end));
eval(['axis([0 ' num2str(xaxis(3)) ' 0 ' num2str(yaxis(3)) '])']);
xlabel('True Strain [-]','FontSize',18);
ylabel('True Stress [MPa]','FontSize',18);
axis square;
box on
end
%-----------------------------------
% Young's Modulus
%-----------------------------------
YoungExp(zz) = getYoungModulus(EngStrain_longi,EngStress,plotname) 
disp(['Young modulus  [GPa] = ' num2str(YoungExp(zz)/1000)]);
figure(8)%figure5=figure('Name','Lankford');
plot(str2num(orientation),YoungExp(zz)/1000,'ok'); hold on
if zz==NbTest
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'XTick',[0:15:90]);
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
% set(gca,'FontSize',18);
% set(gca,'FontSize',18);
xaxis(8)=1.01*max(YoungExp);
yaxis(8)=0.99*min(YoungExp);
set(gca,'YTick',[floor(yaxis(8)/1000):2:ceil(xaxis(8)/1000)]);
eval(['axis([-5 95 ' num2str(yaxis(8)/1000) ' ' num2str(xaxis(8)/1000) '])']);
xlabel('Material Orienation [°]','FontSize',18);
ylabel('Youngs Modulus [GPa]','FontSize',18);
axis square;
box on
end
%-------------------------------------------------------
% Create plastic strain and "plastic" stress variables
%-------------------------------------------------------  
switch YoungId
   case 'Exp'       
        PlasticStrain_longi = TrueStrain_longi - TrueStress(1:length(TrueStrain_longi))./YoungExp(zz);
        PlasticStrain_width = TrueStrain_width + PoissonTheo.*TrueStress(1:length(TrueStrain_width))./YoungExp(zz);
            
   case 'Theo'    
        PlasticStrain_longi = TrueStrain_longi - TrueStress(1:length(TrueStrain_longi))./YoungTheo;
        PlasticStrain_width = TrueStrain_width + PoissonTheo.*TrueStress(1:length(TrueStrain_width))./YoungTheo;
end

%Find the minimum plastic strain by looking at the minimum value from
%the end of the table (avoid local fluctuation)
    minP=length(PlasticStrain_longi);
    for k = minP:-1:1
        if PlasticStrain_longi(k)<=minPStrain
            minP=k+1;
            break
        end
    end

PlasticStrain_longi = PlasticStrain_longi(minP:end);
PlasticStress = TrueStress(minP:end);
PlasticStrain_width = PlasticStrain_width(minP:end);

figure(4)%figure4=figure('Name','True Stress/ Plastic Strain');
plot(PlasticStrain_longi,PlasticStress); hold on
pos=length(PlasticStrain_longi);%randi([length(TrueStrain_longi)-50,length(TrueStrain_longi)],1);
eval(['text(PlasticStrain_longi(pos),PlasticStress(pos),''\leftarrow ' plotname ' '')']);
if zz==zz(end)
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
% set(gca,'FontSize',18);
% set(gca,'FontSize',18);
xaxis(4)=max(xaxis(4),1.1*PlasticStrain_longi(end));
yaxis(4)=max(yaxis(4),1.1*PlasticStress(end));
eval(['axis([0 ' num2str(xaxis(4)) ' 0 ' num2str(yaxis(4)) '])']);
xlabel('Plastic Strain [-]','FontSize',18);
ylabel('True Stress [MPa]','FontSize',18);
axis square;
box on
end
%-------------------------------------------------------
% Calculate Lankford Coefficients 
%------------------------------------------------------- 
lossofcorrelationshort=min(length(PlasticStrain_longi),length(PlasticStrain_width));
% PlasticStrain_longi=PlasticStrain_longi(1:lossofcorrelationshort)
% PlasticStrain_width=PlasticStrain_width(1:lossofcorrelationshort)
% 
% [kk, ia, ic]  = unique(PlasticStrain_longi);
% PlasticStrain_longi = PlasticStrain_longi(ia);
% PlasticStrain_width = PlasticStrain_width(ia);
% 
% [kk, ia, ic]  = unique(PlasticStrain_width);
% PlasticStrain_longi = PlasticStrain_longi(ia);
% PlasticStrain_width = PlasticStrain_width(ia);

PlasticStrain_thick =-PlasticStrain_longi(1:lossofcorrelationshort)-PlasticStrain_width(1:lossofcorrelationshort);
Lankford = fit(PlasticStrain_thick,PlasticStrain_width(1:lossofcorrelationshort),'poly1');
%-------------------------------------------------------
% REMOVE OUTLIERS FOR LANKFORD CALCULATION
%-------------------------------------------------------
fdata = feval(Lankford,PlasticStrain_thick);
I = abs(fdata - PlasticStrain_width(1:lossofcorrelationshort)) > 0.5*std(PlasticStrain_width(1:lossofcorrelationshort)); 
outliers = excludedata(PlasticStrain_thick,PlasticStrain_width(1:lossofcorrelationshort),'indices',I);
Lankford = fit(smooth(PlasticStrain_thick),smooth(PlasticStrain_width(1:lossofcorrelationshort)),'poly1','Exclude',outliers);
coeffnames(Lankford);
clear I fdata outliers 
%-------------------------------------------------------
figure(5)%figure5=figure('Name','Lankford');
plot(PlasticStrain_thick,PlasticStrain_width); hold on
plot(PlasticStrain_thick,Lankford.p1.*PlasticStrain_thick+Lankford.p2); hold on
pos=randi([length(PlasticStrain_thick)-50,length(PlasticStrain_thick)],1);
eval(['text(PlasticStrain_thick(pos),PlasticStrain_width(pos),''\leftarrow ' plotname ' '')']);
% if zz==zz(end)
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
% set(gca,'FontSize',18);
% set(gca,'FontSize',18);
xaxis(5)=min(xaxis(5),1.1*PlasticStrain_thick(end));
yaxis(5)=min(yaxis(5),1.1*PlasticStrain_width(end));
eval(['axis([' num2str(xaxis(5)) ' 0 ' num2str(yaxis(5)) ' 0 ' '])']);

xlabel('Plastic thick. strain  [-]','FontSize',18);
ylabel('Plastic width strain [-]','FontSize',18);
axis square;
box on
% end
%-------------------------------------------------------
% Strain Rate
%------------------------------------------------------- 
TrueStrainRate=diff(TrueStrain_longi)./diff(TrueTime);

figure(6)%figure5=figure('Name','Lankford');
plot(TrueStrain_longi(2:end),smooth(TrueStrainRate)); hold on
pos=randi([length(TrueStrainRate)-10,length(TrueStrainRate)],1);
eval(['text(TrueStrain_longi(pos),TrueStrainRate(pos),''\leftarrow ' plotname ' '')']);
if zz==zz(end)
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
% set(gca,'FontSize',18);
% set(gca,'FontSize',18);
xaxis(6)=max(xaxis(6),1.1*TrueStrain_longi(end));
yaxis(6)=max(yaxis(6),2.0*mean(TrueStrainRate));
eval(['axis([0 ' num2str(xaxis(6)) ' 0 ' num2str(yaxis(6)+yaxis(6)) '])']);
xlabel('True Strain [-]','FontSize',18);
ylabel('True Strain Rate [1/s]','FontSize',18);
axis square;
box on
end

%-------------------------------------------------------
% Save data Experimental Results
%------------------------------------------------------- 
A=zeros(length(EngStrain_longi),12);
A(:,1)=Time(1:length(EngStrain_longi));
A(:,2)=u1(1:length(EngStrain_longi));
A(:,3)=Force(1:length(EngStrain_longi));
A(:,4)=EngStrain_longi(1:length(EngStrain_longi));
A(:,5)=EngStress(1:length(EngStrain_longi));
A(1:length(TrueStrain_longi),6)=TrueStrain_longi;
A(1:length(TrueStress),7)=TrueStress;
A(2:length(TrueStrainRate)+1,8)=TrueStrainRate;
A(1:length(PlasticStrain_longi),9)=PlasticStrain_longi;
A(1:length(PlasticStress),10)=PlasticStress;
A(1:length(PlasticStrain_width),11)=PlasticStrain_thick;
A(1:length(PlasticStrain_thick),12)=PlasticStrain_width;
A(1,13)=Lankford.p1;
A(1,14)=thickness;
A(1,15)=outerwidth;
% A(1,15)=innerwidth;
A(1,16)=YoungExp(zz);

mkdir _results_full
cd _results_full
name = strcat(testname,'.csv');
fid=fopen(name,'w');
fprintf(fid,'Time[s],Displacement[mm],Force[kN],EngStrain_longi,EngStress[MPa],TrueStrain,TrueStress[MPa],TrueStrainRate[1/s],PlasticStrain_longi,PlasticStress[MPa],PlasticStrain_thick,PlasticStrain_width,R-Value,Thickness[mm],Width[mm], Young''s Modulus [MPa]\n');  
fclose(fid);
dlmwrite(name,A,'-append','precision','%.8f');
cd ..
clear A
%-------------------------------------------------------
% Clear Variables 
%------------------------------------------------------- 
clear Time u* Force Eng* True* Plastic* extenso* loss* minP Lankford testname fmax i k ind
clear filename geometry orientation testnumber testname outerwidth innerwidth thickness NominalSpeed Material mmperpix name
end
%-------------------------------------------------------
% CALCULATE AVERAGE YOUNGS MODULUS CASE EXP
%-------------------------------------------------------
% switch YoungId
%    case 'Exp'
    YoungExpAvg=mean(YoungExp);
    disp(['Average Young modulus  [GPa] = ' num2str(YoungExpAvg/1000)]);
%end
%-------------------------------------------------------
cd _results_full
savefig(1,'_UT_Force_Displacement')
saveas(1,['_UT_Force_Displacement.png'],'png');
savefig(2,'_UT_EngStress_EngStrain');
saveas(2,['_UT_EngStress_EngStrain.png'],'png');
savefig(3,'_UT_TrueStress_TrueStrain');
saveas(3,['_UT_TrueStress_TrueStrain.png'],'png');
savefig(4,'_UT_TrueStress_PlasticStrain');
saveas(4,['_UT_TrueStress_PlasticStrain.png'],'png');
savefig(5,'_UT_R-Values');
saveas(5,['_UT_R-Values.png'],'png');
savefig(6,'_UT_TrueStrainRate_TrueStrain');
saveas(6,['_UT_TrueStrainRate_TrueStrain.png'],'png');
savefig(7,'_UT_Time_Displacement');
saveas(7,['_UT_Time_Displacement.png'],'png');
savefig(8,'_Youngs Modulus');
saveas(8,['_Youngs Modulus.png'],'png');
savefig(9,'_Youngs Modulus_fit');
saveas(9,['_Youngs Modulus_fit.png'],'png');
cd ..
hold off






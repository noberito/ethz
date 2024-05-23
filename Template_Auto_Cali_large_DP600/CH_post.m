function CH_post
clear all; close all; clc;
cd _data
%-----------------------------------
% FIND AND LOAD CH FILES
%-----------------------------------
listOf = dir('*.txt');
ii=1;
for i=1:numel(listOf)
str=listOf(i).name;
if (strfind(str,'CH_'))
filenames{ii}=listOf(i).name;
ii=ii+1;
end
end
NbTest = numel(filenames);
clear i ii listOf str
%
xaxis=zeros(1,10);
yaxisax1=zeros(1,10);
yaxisax2=zeros(1,10);
yaxis=zeros(1,10);
cd ..
%-----------------------------------
%Create F/U log strain plot
%-----------------------------------
fig=figure(1);
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
ax1=gca;
set(ax1,'Color','none');
set(ax1,'XMinorTick','on');
% set(ax1,'XTick',[0:15:90]);
set(ax1,'YMinorTick','on');
set(ax1,'TickLength',[0.01 0.01]);
set(ax1,'FontSize',18);
set(ax1,'LineWidth',1);
xlabel('Displacement [mm]','FontSize',18);
ylabel('Force [kN]','FontSize',18);
box(ax1,'off')
axis square
%-----------------------------------
%Create second axis
%-----------------------------------
ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','b');
linkaxes([ax1 ax2], 'x'); 
hold on
axis square
set(ax2,'XTickLabel',[]);
set(ax2,'XMinorTick','on');
% set(ax2,'XTick',[0:15:90],'fontsize',15);
% set(ax2,'YTick',[0:0.05:2.0],'fontsize',15);
set(ax2,'YMinorTick','on');
set(ax2,'TickLength',[0.01 0.01]);
set(ax2,'FontSize',18);
set(ax2,'LineWidth',1);
ylabel('Local Ax. Strain [-]','fontsize',18,'color','b');
%-----------------------------------
% Identify specimen type 
%-----------------------------------
for zz=1:NbTest
cd _data
%-----------------------------------
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
diameter =str2num(cell2mat(A{1}(9))); 
mmperpix=str2num(cell2mat(A{1}(10)));
extenso_zone_long_u1=str2num(cell2mat(A{1}(11)));
extenso_zone_short_u2=str2num(cell2mat(A{1}(12)));
extenso_zone_short_u3=str2num(cell2mat(A{1}(13)));
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
for i=1:3
A(any(isnan(A(:,1)),2),:) = []
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
Time = A(:,1);
u1=A(:,2);
u2=A(:,3);
u3=A(:,4);
% Force = A(:,5);
Force = smooth(A(:,5),'lowess');
clear A;
[fmax ind] = max(Force);
%-----------------------------------
%Plotting Force/Displacement
%-----------------------------------
xaxis(1)=max(xaxis(1),max(u1));
yaxisax1(1)=max(yaxisax1(1),max(Force));
set(1,'CurrentAxes',ax1);
plot(u1,Force,'LineStyle','none','Marker','.','MarkerSize',3,'MarkerEdgeColor','k'); hold on; 
pos=randi([ind-50,ind],1);
eval(['text(u1(pos),Force(pos),''\leftarrow ' plotname ' '')']);
%-----------------------------------
%Plotting Time Displacement
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
plot(Time,u1,'.k','MarkerSize',4); hold on
% pos=randi([length(Time)-200,length(Time)],1);
eval(['text(Time(pos),u1(pos),''\leftarrow ' plotname ' '')']);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
% set(gca,'FontSize',18);
% set(gca,'FontSize',18);
xaxis(2)=max(xaxis(2),1.1*Time(end));
yaxis(2)=max(yaxis(2),1.1*max(u1));
eval(['axis([0 ' num2str(xaxis(2)) ' 0 ' num2str(yaxis(2)) '])']);
xlabel('Time [s]','FontSize',18);
ylabel('Displacement [mm]','FontSize',18);
axis square;
box on
%-----------------------------------
% Find fracture
%-----------------------------------
[fmax ind]=max(Force);
AForce=diff(Force(ind:end));
AU=diff(u1(ind:end));
ATime=diff(Time(ind:end));
AForceATime=AForce./ATime;
AUATime=AU./ATime;

a=-AForceATime;
% ind2 = find(a>0.5);
ind2 = find(a>0.00004);
ind3 = find(abs(AUATime)>0.1);
% if isempty(ind2)==1
%     ind2=1;
% end
% if isempty(ind3)==1
%     ind3=ind2;
% end
% 
% FractureIndForce=ind+ind2(1)-1;
% FractureIndDisp=ind+ind3(1)-1;
% FractureInd = min(FractureIndForce,FractureIndDisp);

% FractureInd=min(FractureInd,min(lossofcorr(1)));
FractureInd=min(length(Force));
%-----------------------------------
figure
plot(u1,Force,'.k'); hold on; 
xlabel('Displacement [mm]','FontSize',18);
ylabel('Force [kN]','FontSize',18);
% eval(['axis([0 ' num2str(xaxis(1)) ' 0 ' num2str(yaxisax1(1)) '])']);
yyaxis right
plot(u1(ind+1:end),AForceATime,'.b');
plot(u1(ind+1:end),AUATime,'.r');
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
ylabel('\DeltaF/\Deltat [kN/s]','FontSize',18);
axis square;
box on
%-----------------------------------
figure(1)
set(1,'CurrentAxes',ax1);
plot(u1(1:FractureInd),Force(1:FractureInd),'r'); hold on; 
figure(2)
plot(Time(1:FractureInd),u1(1:FractureInd),'r'); hold on; 

%-----------------------------------
% Local Strain
%-----------------------------------
LocStrain_longi(:,1)=log(u2./extenso_zone_short_u2+1.);
LocStrain_longi(:,2)=log(u3./extenso_zone_short_u3+1.);  
figure(1);
set(1,'CurrentAxes',ax2);
for i=1:size(LocStrain_longi,2)
yaxisax2(1)=real(max(yaxisax2(1),max(LocStrain_longi(:,i))));
plot(u1,LocStrain_longi(:,i),'.k','MarkerSize',4); hold on; 
plot(u1(1:FractureInd),LocStrain_longi(1:FractureInd,i),'b'); 
% pos=randi([length(u1)-200,length(u1)],1);
eval(['text(u1(pos),LocStrain_longi(pos,i),''\leftarrow ' plotname ' '')']);
end
%-------------------------------------------------------
% Save data
%------------------------------------------------------- 
% A=zeros(length(EngStrain_longi),4);
A(:,1)=Time(1:FractureInd);
A(:,2)=u1(1:FractureInd);
A(:,3)=Force(1:FractureInd);
A(:,4)=LocStrain_longi(1:FractureInd,1);
A(:,5)=LocStrain_longi(1:FractureInd,2);
A(1,6)=thickness;
A(1,7)=outerwidth;
A(1,8)=diameter;

mkdir _results
cd _results
name = strcat(testname,'.csv')
fid=fopen(name,'w');
fprintf(fid,'Time[s],Displacement[mm],Force[kN],AxStrain_1,AxStrain_2,Thickness[mm],OuterWidth[mm],Diameter[mm]\n');  
fclose(fid);
dlmwrite(name,A,'-append');
cd ..
clear A
%-------------------------------------------------------
% Clear Variables 
%------------------------------------------------------- 
clear a A* ind* fmax Force Fract* Loc* extenso* ...
    force u* Mat* mmperpix NominalSpeed thickness outerwidth innerwidth ...
    Time orientation lossofcorr i testnumber filename testname plotname pos zz
end

figure(1)
set(1,'CurrentAxes',ax1);
eval(['axis([0 ' num2str(1.1*xaxis(1)) ' 0 ' num2str(1.1*yaxisax1(1)) '])']);
set(1,'CurrentAxes',ax2);
eval(['axis([0 ' num2str(1.1*xaxis(1)) ' 0 ' num2str(1.3*yaxisax2(1)) '])']);
set(1,'CurrentAxes',ax1);


%-------------------------------------------------------
% Save PLOTS
%------------------------------------------------------- 
cd _results
savefig(1,'_CH_Force_Disp_LocAxStrain_all')
saveas(1,['_CH_Force_Disp_LocAxStrain_all.png'],'png');
savefig(2,'_CH_Time_Displacement_all');
saveas(2,['_CH_Time_Displacement_all.png'],'png');
% savefig(3,'_CH_LocAxStrain_Displacement');
% saveas(3,['_CH_LocAxStrain_Displacement.png'],'png');
cd ..
hold off
%%
%-------------------------------------------------------
% PAUSE
%------------------------------------------------------- 
pause()
%-------------------------------------------------------
% FOR SIMULATION
%-------------------------------------------------------
cd _results
clear all; close all; clc;
listOf = dir('*.csv');
ii=1;
for i=1:numel(listOf)
if (strncmp(listOf(i).name,'CH_',3)==1)
filenames{ii}=listOf(i).name;
ii=ii+1;
end
end
clear i ii listOf
%-----------------------------------
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
a=A(2:end,2);
aa=any(a~=0,2);
A=A(aa,:);

eval(['Time(testnumber).' nameori '=A(:,1);']); 
eval(['u(testnumber).' nameori '=A(:,2);']);   
eval(['Force(testnumber).' nameori '=A(:,3);']);  
eval(['LocStrain_longi1(testnumber).' nameori '=A(:,4);']); 
eval(['LocStrain_longi2(testnumber).' nameori '=A(:,5);']); 
eval(['thickness(testnumber).' nameori '=A(1,6);']);  
eval(['widthOuter(testnumber).' nameori '=A(1,7);']);  
eval(['diameter(testnumber).' nameori '=A(1,8);']);  
clear A a aa name geometry orientation testnumber plotname testname nameori
end
%-------------------------------------------------------
% Normalize Force
%------------------------------------------------------- 
names=fieldnames(Force);
for i=1:numel(names);
xaxis=0;
yaxis=0;
name=char(names(i));   
geometry= name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];
figure(i)
for testnumber=1:100
plotname= [geometry '-' orientation '-' num2str(testnumber)];
try
eval(['area(testnumber).' nameori '=(widthOuter(testnumber).' nameori '-diameter(testnumber).' nameori ').*thickness(testnumber).' nameori ';']);
eval(['Force_area(testnumber).' nameori '=Force(testnumber).' nameori './area(testnumber).' nameori ';']);
eval(['plot(u(testnumber).' nameori ',Force_area(testnumber).' nameori ');']); hold on
eval(['pos=randi([length(u(testnumber).' nameori '-20),length(u(testnumber).' nameori ')],1);']);
eval(['text(u(testnumber).' nameori '(pos),Force_area(testnumber).' nameori '(pos),''\leftarrow ' plotname ' '')']);
eval(['xaxis=max(xaxis,max(abs(u(testnumber).' nameori ')));']);
eval(['yaxis=max(yaxis,max(abs(Force_area(testnumber).' nameori ')));']);
%
eval(['maxDisp(testnumber)=max(u(testnumber).' nameori ');']); 
end
end
%-------------------------------------------------------
eval(['u_int.' nameori '=linspace(0,min(maxDisp),1000);']);
% figure
for testnumber=1:100
    try
    clear temp_F temp_u kk ia
    eval(['[temp_u kk ia]=unique(u(testnumber).' nameori ');']);
    eval(['temp_F=Force(testnumber).' nameori ';']);
    temp_F=temp_F(kk);
    eval(['Force_int(testnumber).' nameori '=pchip(temp_u,temp_F,u_int.' nameori ');']); 
%     eval(['plot(u_int.' nameori ',Force_int(testnumber).' nameori ',''-k'')']); hold on
    if abs(temp_F(1))>0
    ll=testnumber;
    end
    clear temp_F temp_u kk ia
    end
end
eval(['tempF=squeeze(cat(ll,Force_int.' nameori '));']);
eval(['stdd.' nameori '=mean(std(tempF,0,2));']);
% DETERMINE ERROR
%for ii=1:3 %Hardcoded because numbering of tests is off 
for ii=1:ll  
    eval(['err(ii)= mean(abs(tempF(:,ii) - mean(tempF,2)));']);    
    if err(ii) > 2.*mean(std(tempF,0,2))
       usevector(1,ii)=0; 
    else
       usevector(1,ii)=1;
    end         
end 
eval(['[kk best.' nameori ']=min(err);']);
clear kk tempF stdd ii ll maxDisp err
%-------------------------------------------------------
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
% set(gca,'FontSize',18);
% set(gca,'FontSize',18);
xaxis=1.1*xaxis;
yaxis=1.1*yaxis;
eval(['axis([0 ' num2str(xaxis) ' 0 ' num2str(yaxis) '])']);
xlabel('Displacement [mm]','FontSize',18);
ylabel('Force/Area [GPa]','FontSize',18);
axis square;
box on
eval(['savefig(i,''_' nameori '_Force_Disp_normalized'')']);
saveas(1,['CH_Force_Disp_normalized.png'],'png')
clear a aa name geometry orientation testnumber plotname testname nameori A
end
cd ..
%-------------------------------------------------------
% Save data Simulation
%------------------------------------------------------- 
clear A
names=fieldnames(Force);
%testnumber=1
for i=1:numel(names);
name=char(names(i));   
geometry= name(1:2);
orientation= name(4:5);
nameori = [geometry '_' orientation];

eval(['A(:,1)=Time(best.' nameori ').' nameori ';']);
eval(['A(:,2)=u(best.' nameori ').' nameori ';']);
eval(['A(:,3)=Force(best.' nameori ').' nameori ';']);
eval(['A(:,4)=LocStrain_longi1(best.' nameori ').' nameori ';']);
eval(['A(:,5)=LocStrain_longi2(best.' nameori ').' nameori ';']);
eval(['A(1,6)=thickness(best.' nameori ').' nameori ';']);
eval(['A(1,7)=widthOuter(best.' nameori ').' nameori ';']);
eval(['A(1,8)=diameter(best.' nameori ').' nameori ';']);
eval(['A(1,9)=best.' nameori ';']);

try
    mkdir _simulation
end
cd _simulation
name = strcat(nameori,'_exp.csv')
dlmwrite(name,A);
cd ..
disp('CHOSEN FOR SIMULATION:')
disp(best)
clear A
end


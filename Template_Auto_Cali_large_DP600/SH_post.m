function SH_post
clear all; close all; clc;
cd _data
%-----------------------------------
% FIND AND LOAD SH FILES
%-----------------------------------
listOf = dir('*.txt');
ii=1;
for i=1:numel(listOf)
str=listOf(i).name;
if (strfind(str,'SH_'))
filenames{ii}=listOf(i).name;
ii=ii+1;
end
end
NbTest = numel(filenames);
clear i ii listOf str
%
xaxis=zeros(1,10);
yaxis=zeros(1,10);
cd ..
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

if (strcmp(orientation,'-45'))
    orientation='m45'
elseif (strcmp(orientation,'45'))  
    orientation='p45'
else
    orientation=strcat('0', orientation);
end
   
testnumber=cell2mat(A{1}(5));
NominalSpeed=str2num(cell2mat(A{1}(6)));
thickness =str2num(cell2mat(A{1}(7)));%2.55;
widthOuter =str2num(cell2mat(A{1}(8))); %10;
widthInner =str2num(cell2mat(A{1}(9))); 
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
for zz=1:3
A(any(isnan(A(:,1)),2),:) = [];
a=A(:,zz+1);
aa=any(a~=-1,2);
B=a(aa,:);
lossofcorr(zz)=length(B);
clear a aa B zz
end   
%-----------------------------------
% Force displacement 
%-----------------------------------
clear Time u1 u2 u3 EffStrain
Time = A(:,1);
u1=A(:,2);
u2=A(:,3);
u3=A(:,4);
%EffStrain=A(:,19);

count = 1;
for iii=7:2:10
    ExtensoX(:,count)=A(:,iii);
    ExtensoY(:,count)=A(:,iii+1);
    count=count+1;
end
clear count iii
count = 1;
for iii=1:2:size(ExtensoX,2)
Ux(:,count) = abs( (ExtensoX(:,iii+1)-ExtensoX(:,iii))-(ExtensoX(1,iii+1)-ExtensoX(1,iii)) ).*mmperpix;
Uy(:,count) = abs( (ExtensoY(:,iii+1)-ExtensoY(:,iii))-(ExtensoY(1,iii+1)-ExtensoY(1,iii)) ).*mmperpix;
count=count+1;
end 
clear count iii


% Force = A(:,5);
Force = smooth(A(:,5),'lowess');
clear A;
[fmax ind] = max(Force);
%---------------------------------
figure(1);
plot(u1,Force,'.k','MarkerSize',4); hold on; 
pos=randi([ind-5,ind],1);
eval(['text(u1(pos),Force(pos),''\leftarrow ' plotname ' '')']);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18);
set(gca,'LineWidth',1);
% set(gca,'FontSize',18);
% set(gca,'FontSize',18);
xaxis(1)=max(xaxis(1),1.1*max(u1));
yaxis(1)=max(yaxis(1),1.1*max(Force));
eval(['axis([0 ' num2str(xaxis(1)) ' 0 ' num2str(yaxis(1)) '])']);
xlabel('Displacement [mm]','FontSize',18);
ylabel('Force [kN]','FontSize',18);
axis square;
box on
%-----------------------------------
figure(2);
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
ind2 = min(find(a>0.0005));
% ind2 = find(a>0.001);

% ind3 = min(find(AUATime<0.0)); 
ind3 = min(find(abs(AUATime)>0.01));

if isempty(ind3)==1
    ind3=ind2;
end

if isempty(ind2)==1
    ind2=0;
end
if isempty(ind3)==1
    ind3=0;
end

FractureIndForce=ind+ind2(1)-1;
FractureIndDisp=min(ind+ind3(1)-1);
FractureInd = min(FractureIndForce,FractureIndDisp);

FractureInd=min(FractureInd,min(lossofcorr(1)));
%-----------------------------------
% figure
% plot(u1,Force,'-k'); hold on; 
% xlabel('Displacement [mm]','FontSize',18);
% ylabel('Force [kN]','FontSize',18);
% eval(['axis([0 ' num2str(xaxis(1)) ' 0 ' num2str(yaxis(1)) '])']);
% yyaxis right
% plot(u1(ind+1:end),AForceATime,'.b');
% set(gca,'XMinorTick','on');
% set(gca,'YMinorTick','on');
% set(gca,'TickLength',[0.01 0.01]);
% set(gca,'FontSize',18);
% set(gca,'LineWidth',1);
% ylabel('\DeltaF/\Deltat [kN/s]','FontSize',18);
% axis square;
% box on
%-----------------------------------
figure(1)
plot(u1(1:FractureInd),Force(1:FractureInd),'r'); hold on; 
figure(2)
plot(Time(1:FractureInd),u1(1:FractureInd),'r'); hold on; 
%-------------------------------------------------------
% Save data
%------------------------------------------------------- 
A(:,1)=Time(1:FractureInd);
A(:,2)=u1(1:FractureInd);
A(:,3)=Force(1:FractureInd);
%A(:,4)=EffStrain(1:FractureInd);%Uy(1:FractureInd);
A(:,4)=0.;%Uy(1:FractureInd);
A(:,5)=0.;
A(1,6)=thickness;
A(1,7)=widthOuter;
A(1,8)=widthInner;

mkdir _results
cd _results
name = strcat(testname,'.csv')
fid=fopen(name,'w');
fprintf(fid,'Time[s],Displacement longi[mm],Force[kN],Effective Strain,N/A,Thickness[mm],OuterWidth[mm],N/A\n');  %,Displacement width[mm]
fclose(fid);
dlmwrite(name,A,'-append');
cd ..
clear A
%-------------------------------------------------------
% Clear Variables 
%------------------------------------------------------- 
clear a A* ind* fmax Force Fract* Loc* extenso* ...
    force u* Mat* mmperpix Speed thickness width* Time Extenso* U*
end
%-------------------------------------------------------
% Save PLOTS
%------------------------------------------------------- 
cd _results
savefig(1,'_SH_Force_Displacement')
saveas(1,['_SH_Force_Displacement.png'],'png');
savefig(2,'_SH_Time_Displacement');
saveas(2,['_SH_Time_Displacement.png'],'png');
cd ..
% clear all
hold off
%
%-------------------------------------------------------
%PAUSE
%------------------------------------------------------- 
pause()
%-------------------------------------------------------
% EffStrainFiles
%------------------------------------------------------- 
% cd _results
% clear all; close all; clc;
% listOf = dir('*.csv');
% ii=1;
% for i=1:numel(listOf)
% if (strncmp(listOf(i).name,'SH',2)==1)
% filenames{ii}=listOf(i).name;
% ii=ii+1;
% end
% end
% %clear i ii listOf
% %-----------------------------------
% for i=1:numel(filenames)
% name = char(filenames(i))
% geometry =  name(1:4);
% orientation = name(6:8);
% testnumber  = name(10);
% plotname = [geometry '-' orientation '-' testnumber];
% plotnames{i} = plotname;
% testname = [geometry '_' orientation '_' testnumber];    
% nameori = [geometry '_' orientation];
% testnumber = str2num(name(10));
% 
% A=csvread(char(filenames(i)),1,0);
% a=A(2:end,2);
% aa=any(a~=0,2);
% A=A(aa,:);
% 
% Time(:,i)=[A(:,1);zeros(5000-size(A,1),1)];
% u(:,i)=[A(:,2);zeros(5000-size(A,1),1)];
% Force(:,i)=[A(:,3);zeros(5000-size(A,1),1)];
% EffStrain(:,i)=[A(:,4);zeros(5000-size(A,1),1)];
% 
% clear A a aa name geometry orientation testnumber plotname testname nameori
% end
% %-----------------------------------
% count=1
% %-----------------------------------
% % 00° Direction
% %-----------------------------------
% for i=1:numel(filenames)
% str=plotnames{i};
% if strfind(str,'SHLD')
%     if strfind(str,'000')
%         x1(:,count) = EffStrain(:,i);
%         names{count}=str;
%         count=count+1;
%     end 
% end
% end
% 
% for i=1:numel(filenames)
% str=plotnames{i};
% if strfind(str,'SHMD')
%     if strfind(str,'000')
%         x1(:,count) = EffStrain(:,i);
%         names{count}=str;
%         count=count+1;
%     end 
% end
% end
% 
% for i=1:numel(filenames)
% str=plotnames{i};
% if strfind(str,'SHHD')
%     if strfind(str,'000')
%         x1(:,count) = EffStrain(:,i);
%         names{count}=str;
%         count=count+1;
%     end 
% end
% end
% %-----------------------------------
% % 45° Direction
% %-----------------------------------
% for i=1:numel(filenames)
% str=plotnames{i};
% if strfind(str,'SHLD')
%     if strfind(str,'45')
%         x1(:,count) = EffStrain(:,i);
%         names{count}=str;
%         count=count+1;
%     end 
% end
% end
% 
% for i=1:numel(filenames)
% str=plotnames{i};
% if strfind(str,'SHMD')
%     if strfind(str,'45')
%         x1(:,count) = EffStrain(:,i);
%         names{count}=str;
%         count=count+1;
%     end 
% end
% end
% 
% for i=1:numel(filenames)
% str=plotnames{i};
% if strfind(str,'SHHD')
%     if strfind(str,'45')
%         x1(:,count) = EffStrain(:,i);
%         names{count}=str;
%         count=count+1;
%     end 
% end
% end
% %-----------------------------------
% % 90° Direction
% %-----------------------------------
% for i=1:numel(filenames)
% str=plotnames{i};
% if strfind(str,'SHLD')
%     if strfind(str,'090')
%         x1(:,count) = EffStrain(:,i);
%         names{count}=str;
%         count=count+1;
%     end 
% end
% end
% 
% for i=1:numel(filenames)
% str=plotnames{i};
% if strfind(str,'SHMD')
%     if strfind(str,'090')
%         x1(:,count) = EffStrain(:,i);
%         names{count}=str;
%         count=count+1;
%     end 
% end
% end
% 
% for i=1:numel(filenames)
% str=plotnames{i};
% if strfind(str,'SHHD')
%     if strfind(str,'090')
%         x1(:,count) = EffStrain(:,i);
%         names{count}=str;
%         count=count+1;
%     end 
% end
% end
% clear count
% %-----------------------------------
% fig=figure(1);
% set(fig,'color','w'); 
% set(fig,'PaperSize', [5 5]);
% hold on
% x0=0;
% y0=0;
% width=600;
% height=660;
% set(fig,'units','points','position',[0,0,width,height])
% 
% for i=1:numel(filenames)
%     temp = x1(:,i);     
%     aa=any(temp~=0,2);
%     temp=temp(aa,:);
%     color = jet(size(temp,1));
%     for ii=1:length(temp)
%         plot(i,temp(ii),'ok', 'MarkerSize',5, 'MarkerEdgeColor', color(ii,:), 'MarkerFaceColor',  color(ii,:)); hold on
%     end
%     clear temp color aa ii
% end
% 
% yaxismin=min(min(x1(x1>0)));
% yaxismax=max(max(x1));
% 
% 
% ax1=gca;
% set(ax1,'Color','none');
% set(ax1,'XMinorTick','off');
% set(ax1,'XTick',[0:1:18]);
% set(ax1,'YMinorTick','on');
% set(ax1,'TickLength',[0.01 0.01]);
% set(ax1,'FontSize',18);
% set(ax1,'LineWidth',1);
% ylabel('Effective Strain [-]','FontSize',18);
% set(ax1,'XTick',1:18,'XTickLabel',names)
% xtickangle(90)
% axis square;
% box on
% 
% 
% savefig(1,'_SH_Eff_Strain')
% saveas(1,['_SH_Eff_Strain.png'],'png');
% 
% hold off
% cd ..
%-------------------------------------------------------
%PAUSE
%------------------------------------------------------- 
%pause()
%-------------------------------------------------------
% FOR SIMULATION
%------------------------------------------------------- 
cd _results
clear all; close all; clc;
listOf = dir('*.csv');
ii=1;
for i=1:numel(listOf)
if (strncmp(listOf(i).name,'SHMD_',3)==1)
filenames{ii}=listOf(i).name;
ii=ii+1;
end
end
clear i ii listOf
%-----------------------------------
for i=1:numel(filenames)
name=char(filenames(i))
geometry=  name(1:4);
orientation= name(6:8);
testnumber = name(10);
plotname = [geometry '-' orientation '-' testnumber];
testname = [geometry '_' orientation '_' testnumber];    
nameori = [geometry '_' orientation];
testnumber = str2num(name(10));

A=csvread(char(filenames(i)),1,0);
a=A(2:end,2);
aa=any(a~=0,2);
A=A(aa,:);

eval(['Time(testnumber).' nameori '=A(:,1);']); 
eval(['u(testnumber).' nameori '=A(:,2);']);   
eval(['Force(testnumber).' nameori '=A(:,3);']);
eval(['uy(testnumber).' nameori '=A(:,4);']);   
% eval(['LocStrain_longi1(testnumber).' nameori '=A(:,4);']); 
% eval(['LocStrain_longi2(testnumber).' nameori '=A(:,5);']); 
eval(['thickness(testnumber).' nameori '=A(1,6);']);  
eval(['widthOuter(testnumber).' nameori '=A(1,7);']);  
clear A a aa name geometry orientation testnumber plotname testname nameori
end
%-------------------------------------------------------
% NORMALIZE FORCE / REMOVE OUTLIERS / CHOOSE EXP FOR SIM
%------------------------------------------------------- 
names=fieldnames(Force);
for i=1:numel(names);
xaxis=0;
yaxis=0;
name=char(names(i));   
geometry= name(1:4);
orientation= name(6:8);
nameori = [geometry '_' orientation];
figure(i)
for testnumber=1:100
plotname= [geometry '-' orientation '-' num2str(testnumber)];
try
eval(['area(testnumber).' nameori '=widthOuter(testnumber).' nameori '.*thickness(testnumber).' nameori ';']);
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
clear a aa name geometry orientation testnumber plotname testname nameori A
end
cd ..
%-------------------------------------------------------
% Save data Simulation
%------------------------------------------------------- 
clear A
names=fieldnames(Force);
for i=1:numel(names);
name=char(names(i));   
geometry= name(1:4);
orientation= name(6:8);
nameori = [geometry '_' orientation];

eval(['A(:,1)=Time(best.' nameori ').' nameori ';']);
eval(['A(:,2)=u(best.' nameori ').' nameori ';']);
eval(['A(:,3)=Force(best.' nameori ').' nameori ';']);
% eval(['A(:,4)=uy(best.' nameori ').' nameori ';']);
% eval(['A(:,4)=LocStrain_longi1(best.' nameori ').' nameori ';']);
% eval(['A(:,5)=LocStrain_longi2(testnumber).' nameori ';']);
% A(1:length(A(:,4)),5)=0.0;
eval(['A(1,6)=thickness(best.' nameori ').' nameori ';']);
eval(['A(1,7)=widthOuter(best.' nameori ').' nameori ';']);
% eval(['A(1,8)=widthInner(best.' nameori ').' nameori ';']);
eval(['A(1,9)=best.' nameori ';']);

try
    mkdir _simulation
end
cd _simulation
newStr= [geometry(1:2) '_' orientation];
name = strcat(newStr,'_exp.csv')
dlmwrite(name,A);
cd ..
disp('CHOSEN FOR SIMULATION:')
disp(best)
clear A
end
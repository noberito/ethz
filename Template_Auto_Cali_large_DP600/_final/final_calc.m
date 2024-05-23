clear all; close all; format long; clc
global a0 P Pi Sigma K0;
delete screen_global.txt
diary screen_global.txt
scalfact = 1; % Scaling compared to the bigest set of the 3 specimens (big: 1, medium: 2, small:4)
scalfact = scalfact*5;
%%
fid=fopen('material.dat','r');
AA = textscan(fid,'%s',2,'Delimiter',',');
A = textscan(fid,'%s',20,'Delimiter',',');
fclose(fid);
%-----------------------------------
% LOAD AVAILABLE EXPERIMENTS & CALC OFFSET IN FU
%-----------------------------------
listOf = dir('*.csv');
ii=0;
for i=1:numel(listOf)
files{i}=listOf(i).name;
name=char(files(i));
name=name(1:end-4);
if any(regexp(name,'_exp$'))==1
 eval(['experiments.' name(1:end-4) '=csvread(files{i});']);
 % CALC OFFSET IN FU curve
 if any(regexp(name,'UT_'))~=1
    temp=csvread(files{i});
    tempU=temp(:,2);
    tempF=temp(:,3);
    tempoffset = FUoffset(tempU,tempF);
    eval(['deltaU.' name(1:end-4) ' = tempoffset;']);  
    fid = fopen('_Exp_offset.csv','a');
    fprintf(fid,'%s',[name(1:end-4)]);
    fprintf(fid,',');
    fprintf(fid,'%.6f',[tempoffset]);
    fprintf(fid,'\n');
    fid = fclose(fid);
    clear temp*
 end
 ii=ii+1;
 nbExps = ii;
end
end

clear i ii listOf files name

launchcount = 0;
names=fieldnames(experiments);
%---------------------------------------
% Cluster folder cleanup commands
%---------------------------------------
delete lsf*.*;
delete abaqus*.*;
delete temp*.*;
delete *.odb;
%---------------------------------------
% COPY SIMs to TEMPSIM files % LAUNCH COMPUTATIONS
%---------------------------------------
for i=1:nbExps
name=char(names(i))
%---------------------------------------
% UTs
%---------------------------------------
if (strncmp(name,'UT_',3)==1)
  eval(['copyfile(''' name '.inp'',''temp_' name '.inp'')']);
  %exe_command = ['!bsub -n 1 -W 2:00 -R "ib rusage[mem=500,scratch=500]" "abaqus job=temp_' name ' user=YLD2000_3D_explicit.f double cpus=1 memory="500mb" scratch=\$TMPDIR"'];
  exe_command = ['!sbatch -n 1 -t 0-2 --mem-per-cpu=2G --tmp=2G --wrap "abaqus job=temp_' name ' user=VUMAT_explicit.f double cpus=1 scratch=\$TMPDIR"'];
  eval(exe_command); 
  launchcount = launchcount + 1;
%---------------------------------------
% NT20 & NT6
%---------------------------------------
elseif (strncmp(name,'NT',2)==1)
  eval(['copyfile(''' name '.inp'',''temp_' name '.inp'')']);
  %CREATE PARAS FILE
  eval(['tempdata = experiments. ' name ';']);
  THICKNESS	= tempdata(1,6);
  WIDTH = tempdata(1,8);
  eval(['DISPL=tempdata(end,2)-deltaU.' name '+0.0001;']);
  DTIME=tempdata(end,1);
  TRAMP1= 0.01*DTIME
  RAMP1= 0.001*DISPL
  eval(['fid=fopen(''' name '_paras.inp'',''wt'');']);
  %
  DT	 = DTIME/8.E5;
  fprintf(fid,'*PARAMETER\n');
  fprintf(fid,'DISPL = %f\n',DISPL/2.);
  fprintf(fid,'THICKNESS = %f\n',THICKNESS/2.);
  fprintf(fid,'WIDTH = %f\n',WIDTH/10.);
  fprintf(fid,'DTIME = %f\n',DTIME);
  fprintf(fid,'TRAMP1 = %f\n',TRAMP1);
  fprintf(fid,'RAMP1 = %f\n',RAMP1);
  fprintf(fid,'DT = %f\n',DT);
  fclose(fid); 
  %
  %exe_command = ['!bsub -n 12 -W 2:00 -R "ib rusage[mem=300,scratch=300]" "abaqus job=temp_' name ' user=YLD2000_3D_explicit.f double cpus=12 memory="3600mb" scratch=\$TMPDIR"'];
  exe_command = ['!sbatch -n 12 -t 0-2 --mem-per-cpu=4G --tmp=30G --wrap "abaqus job=temp_' name ' user=VUMAT_explicit.f double cpus=12 memory="3600mb" scratch=\$TMPDIR"'];
  eval(exe_command); 
  launchcount = launchcount + 1;
%---------------------------------------
% CENTRAL HOLE CH
%---------------------------------------
elseif (strncmp(name,'CH',2)==1)
  eval(['copyfile(''' name '.inp'',''temp_' name '.inp'')']);
  %CREATE PARAS FILE
  eval(['tempdata = experiments. ' name ';']);
  THICKNESS	= tempdata(1,6);
  WIDTH = tempdata(1,7)/20.;
  eval(['DISPL=tempdata(end,2)-deltaU.' name '+0.0001;']);
  DTIME=tempdata(end,1);
  TRAMP1= 0.01*DTIME
  RAMP1= 0.001*DISPL
  eval(['fid=fopen(''' name '_paras.inp'',''wt'');']);
  %
  DT	 = DTIME/8.E5;
  fprintf(fid,'*PARAMETER\n');
  fprintf(fid,'DISPL = %f\n',DISPL/2.);
  fprintf(fid,'THICKNESS = %f\n',THICKNESS);
  fprintf(fid,'DTIME = %f\n',DTIME);
  fprintf(fid,'DT = %f\n',DT);
  fprintf(fid,'TRAMP1 = %f\n',TRAMP1);
  fprintf(fid,'RAMP1 = %f\n',RAMP1);
  
  if abs(tempdata(1,8))<=6.5
   copyfile('CH_r2_5.geo','CH.geo','f');
   copyfile('CH_final_r2_5.py','CH_final.py','f');
   fprintf(fid,'WIDTH = %f\n',WIDTH);
  end
  if abs(tempdata(1,8))>6.5
   copyfile('CH_r4.geo','CH.geo','f');
   copyfile('CH_final_r4.py','CH_final.py','f');
   fprintf(fid,'WIDTH = %f\n',8./8.);
  end
  fclose(fid);
  %
  %exe_command = ['!bsub -n 12 -W 2:00 -R "ib rusage[mem=300,scratch=300]" "abaqus job=temp_' name ' user=YLD2000_3D_explicit.f double cpus=12 memory="3600mb" scratch=\$TMPDIR"'];
  exe_command = ['!sbatch -n 12 -t 0-2 --mem-per-cpu=4G --tmp=30G --wrap "abaqus job=temp_' name ' user=VUMAT_explicit.f double cpus=12 memory="3600mb" scratch=\$TMPDIR"'];
  eval(exe_command); 
  launchcount = launchcount + 1;
%---------------------------------------
% SHEAR
%---------------------------------------
elseif (strncmp(name,'SH',2)==1)
  eval(['copyfile(''' name '.inp'',''temp_' name '.inp'')']);
  %CREATE PARAS FILE
  eval(['tempdata = experiments. ' name ';']);
  THICKNESS	= tempdata(1,6);
  WIDTH = tempdata(1,7)/20.;% tempdata(1,8);
  eval(['DISPL=tempdata(end,2)+0.001;']);
  % eval(['DISPL=tempdata(end,2)-deltaU.' name '+0.0001;']);
  DTIME=tempdata(end,1); 
  TRAMP1= 0.01*DTIME
  RAMP1= 0.001*DISPL
  eval(['fid=fopen(''' name '_paras.inp'',''wt'');']);
  %
  DT	 = DTIME/6.E5;
  fprintf(fid,'*PARAMETER\n');
  fprintf(fid,'DISPL = %f\n',DISPL);
  fprintf(fid,'THICKNESS = %f\n',THICKNESS);
  fprintf(fid,'WIDTH = %f\n',WIDTH);
  fprintf(fid,'DTIME = %f\n',DTIME);
  fprintf(fid,'TRAMP1 = %f\n',TRAMP1);
  fprintf(fid,'RAMP1 = %f\n',RAMP1);
  fprintf(fid,'DT = %f\n',DT);
  fclose(fid); 
  %
  %exe_command = ['!bsub -n 12 -W 24:00 -R "ib rusage[mem=1000,scratch=1000]" "abaqus job=temp_' name ' user=YLD2000_3D_explicit.f double cpus=12 memory="12000mb" scratch=\$TMPDIR"'];
  exe_command = ['!sbatch -n 12 -t 0-24 --mem-per-cpu=4G --tmp=30G --wrap "abaqus job=temp_' name ' user=VUMAT_explicit.f double cpus=12 memory="12000mb" scratch=\$TMPDIR"'];
 
  eval(exe_command); 
  launchcount = launchcount + 1;
end
clear name
end
clear i ii listOf THICKNESS WIDTH DISPL DTIME tempdata
eval(['disp('' MESSAGE: ' num2str(launchcount) ' Abaqus simulations submitted'')']);
%---------------------------------------
% WAIT FOR ALL TO BE LAUNCHED
%---------------------------------------
pause(5);
while (numel(dir('*.sta'))~=launchcount);
   pause(5);
end
%---------------------------------------
% WAIT FOR ALL TO BE FINISHED AND POSTPROCESS
%---------------------------------------
while (numel(dir('*.lck'))~=0);
 for i=1:nbExps
 name=char(names(i));  
 namelck=['temp_' char(names(i)) '.lck'];
 nameodb=['temp_' char(names(i)) '.odb'];
 nameout=['temp_' char(names(i)) '.out'];
    if ((2-exist(namelck,'file')) && exist(nameodb,'file') && (2-exist(nameout,'file')))
     eval(['copyfile(''temp_' name '.odb'',''temp.odb'')']);
     eval(['copyfile(''temp_' name '.odb'',''' name '.odb'')']);
     pause(5);
     %---------------------------------------
     if (strncmp(name,'UT_',3)==1)  
     !abaqus cae noGUI=UT.py
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'NT20_',5)==1)  
     !abaqus cae noGUI=NT20_final.py  
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'NT6_',4)==1)  
     !abaqus cae noGUI=NT6_final.py  
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'CH_',3)==1)  
     !abaqus cae noGUI=CH_final.py 
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'SH_',3)==1)  
     !abaqus cae noGUI=SH_final.py 
         
     inp = 'SH_final.geo';
     NODES = getnodes(inp);
     ELEMENTS = getelements(inp);
%      elset = '*Elset, elset=BOUND'; %find elements on boundary of gauge section
%      boundEls = getelementset(elset,inp);
%      elset = '*Elset, elset=REF_AREA'; %find elements in ref_area/central area of gauge section
%      ref_area_Els = getelementset(elset,inp);
     elset = '*Elset, elset=POST_AREA'; %find elements in ref_area/central area of gauge section
     refEls = getelementset(elset,inp);
     centerEls=refEls;
     clear refEls elset;

    % Create output + picture
    % write python script
    delete ShearPost.py
    namezz = 'ShearPost.py';
    fid=fopen(namezz,'w');
    %
    fprintf(fid,'from abaqus import *\n');
    fprintf(fid,'from abaqusConstants import *\n');
    fprintf(fid,'import __main__\n');
    fprintf(fid,'import visualization\n');
    fprintf(fid,'import xyPlot\n');
    fprintf(fid,'import displayGroupOdbToolset as dgo\n');
    fprintf(fid,'o1 = session.openOdb(name=''./temp.odb'', readOnly=False)\n');
    fprintf(fid,'session.viewports[''Viewport: 1''].setValues(displayedObject=o1)\n');
    fprintf(fid,'session.viewports[''Viewport: 1''].view.setValues(session.views[''Front''])\n');

    for jj=1:50 % change to 100   
    fprintf(fid,'session.viewports[''Viewport: 1''].odbDisplay.setFrame(step=0, frame=%i)\n', [jj]);
    fprintf(fid,'odb = session.odbs[''./temp.odb''] \n');
    fprintf(fid,'nf = NumberFormat(numDigits=6, precision=0, format=SCIENTIFIC)\n');
    fprintf(fid,'session.fieldReportOptions.setValues(printTotal=OFF, printMinMax=OFF,\n'); 
    fprintf(fid,'    numberFormat=nf) \n');
    fprintf(fid,'session.writeFieldReport(fileName=''abaqus%i.rpt'', append=OFF,\n', [jj]);
    fprintf(fid,'    sortItem=''Element Label'', odb=odb, step=0, frame=%i,\n', [jj]);
    fprintf(fid,'    outputPosition=INTEGRATION_POINT, variable=((''SDV_TRIAX'',INTEGRATION_POINT),\n');
    fprintf(fid,'    (''SDV_LODE'',INTEGRATION_POINT), (''SDV_EQPS'',INTEGRATION_POINT)), ) \n');
    end

    fclose(fid);
    % run python script
    !abaqus cae noGUI=ShearPost.py
    pause(1);
    clear namezz filename fid jj
    delete ShearPost.py
   
    NewNumInc = 0;
    %Get Data from Output
    clear allkk tempt A output delimiterIn headerlinesIn  
    for kkk=1:200%change to 100
        try
        clear A output delimiterIn headerlinesIn filename;
        eval(['filename = ''abaqus' num2str(kkk) '.rpt'';']);
            if exist(filename,'file')
            NewNumInc = NewNumInc + 1;
            delimiterIn = ' ';
            headerlinesIn = 19;
            A = importdata(filename,delimiterIn,headerlinesIn);
            output = A.data;

        %CALCULATE TRIAX/LODE
        
            for j=1:length(output)
                allkk(j,1,kkk) = output(j,1); % ELEMENT
                allkk(j,2,kkk) = output(j,2); % INT POINT  
                allkk(j,3,kkk) = output(j,3); % TRIAX  
                allkk(j,4,kkk) = output(j,4); % LODE 
                allkk(j,5,kkk) = output(j,5); % PEEQ 
            end
            end
        clear A output delimiterIn headerlinesIn filename;
        end 
    end   
    clear kkk j NewNumInc 
     %------------------------------------------------
    %FIND CRITICAL ELEMENT in PostArea and save

    centerEls = allkk(centerEls(:,1),:,:);
    [kk ia] = max(centerEls(:,5,end));
    ElMax = squeeze(centerEls(ia,:,:))';
    eval(['csvwrite(''' name '_CritElOut.out'',ElMax);']);

    clear allkk ia kk ElMax
    eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
    end
    %---------------------------------------
    eval(['copyfile(''temp.out'',''' name '.out'')']);
    eval(['copyfile(''temp.out'',''temp_' name '.out'')']);
    delete temp.odb temp.out abaqus*   
    end
 end
 clear name namelck nameodb nameout
 end

%---------------------------------------
% CHECK POSTPROCESS - MAKE SURE
%---------------------------------------
for i=1:nbExps
 name=char(names(i));  
 namelck=['temp_' char(names(i)) '.lck'];
 nameodb=['temp_' char(names(i)) '.odb'];
 nameout=['temp_' char(names(i)) '.out'];
    if (2-exist(nameout,'file') && exist(nameodb,'file'))
     eval(['copyfile(''temp_' name '.odb'',''temp.odb'')']);
     eval(['copyfile(''temp_' name '.odb'',''' name '.odb'')']);
     pause(1);
     %---------------------------------------
     if (strncmp(name,'UT_',3)==1)  
     !abaqus cae noGUI=UT.py
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'NT20_',5)==1)  
     !abaqus cae noGUI=NT20_final.py  
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'NT6_',4)==1)  
     !abaqus cae noGUI=NT6_final.py  
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'CH_',3)==1)  
     !abaqus cae noGUI=CH_final.py 
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'SH_',3)==1)  
     delete abaqus*.*  
     !abaqus cae noGUI=SH_final.py
     inp = 'SH_final.geo';
     NODES = getnodes(inp);
     ELEMENTS = getelements(inp);
%      elset = '*Elset, elset=BOUND'; %find elements on boundary of gauge section
%      boundEls = getelementset(elset,inp);
%      elset = '*Elset, elset=REF_AREA'; %find elements in ref_area/central area of gauge section
%      ref_area_Els = getelementset(elset,inp);
     elset = '*Elset, elset=POST_AREA'; %find elements in ref_area/central area of gauge section
     refEls = getelementset(elset,inp);
     centerEls=refEls;
     clear refEls elset;

    % Create output + picture
    % write python script
    delete ShearPost.py
    namezz = 'ShearPost.py';
    fid=fopen(namezz,'w');
    %
    fprintf(fid,'from abaqus import *\n');
    fprintf(fid,'from abaqusConstants import *\n');
    fprintf(fid,'import __main__\n');
    fprintf(fid,'import visualization\n');
    fprintf(fid,'import xyPlot\n');
    fprintf(fid,'import displayGroupOdbToolset as dgo\n');
    fprintf(fid,'o1 = session.openOdb(name=''./temp.odb'', readOnly=False)\n');
    fprintf(fid,'session.viewports[''Viewport: 1''].setValues(displayedObject=o1)\n');
    fprintf(fid,'session.viewports[''Viewport: 1''].view.setValues(session.views[''Front''])\n');

    for jj=1:50 % change to 100   
    fprintf(fid,'session.viewports[''Viewport: 1''].odbDisplay.setFrame(step=0, frame=%i)\n', [jj]);
    fprintf(fid,'odb = session.odbs[''./temp.odb''] \n');
    fprintf(fid,'nf = NumberFormat(numDigits=6, precision=0, format=SCIENTIFIC)\n');
    fprintf(fid,'session.fieldReportOptions.setValues(printTotal=OFF, printMinMax=OFF,\n'); 
    fprintf(fid,'    numberFormat=nf) \n');
    fprintf(fid,'session.writeFieldReport(fileName=''abaqus%i.rpt'', append=OFF,\n', [jj]);
    fprintf(fid,'    sortItem=''Element Label'', odb=odb, step=0, frame=%i,\n', [jj]);
    fprintf(fid,'    outputPosition=INTEGRATION_POINT, variable=((''SDV_TRIAX'',INTEGRATION_POINT),\n');
    fprintf(fid,'    (''SDV_LODE'',INTEGRATION_POINT), (''SDV_EQPS'',INTEGRATION_POINT)), ) \n');
    end

    fclose(fid);
    % run python script
    !abaqus cae noGUI=ShearPost.py
    pause(1);
    clear namezz filename fid jj
    delete ShearPost.py
   
    NewNumInc = 0;
    %Get Data from Output
    clear allkk tempt A output delimiterIn headerlinesIn  
    for kkk=1:200%change to 100
        try
        clear A output delimiterIn headerlinesIn filename;
        eval(['filename = ''abaqus' num2str(kkk) '.rpt'';']);
            if exist(filename,'file')
            NewNumInc = NewNumInc + 1;
            delimiterIn = ' ';
            headerlinesIn = 19;
            A = importdata(filename,delimiterIn,headerlinesIn);
            output = A.data;

        %CALCULATE TRIAX/LODE
        
            for j=1:length(output)
                allkk(j,1,kkk) = output(j,1); % ELEMENT
                allkk(j,2,kkk) = output(j,2); % INT POINT  
                allkk(j,3,kkk) = output(j,3); % TRIAX  
                allkk(j,4,kkk) = output(j,4); % LODE 
                allkk(j,5,kkk) = output(j,5); % PEEQ 
            end
            end
        clear A output delimiterIn headerlinesIn filename;
        end 
    end   
    clear kkk j NewNumInc 
     %------------------------------------------------
    %FIND CRITICAL ELEMENT in PostArea and save

    centerEls = allkk(centerEls(:,1),:,:);
    [kk ia] = max(centerEls(:,5,end));
    ElMax = squeeze(centerEls(ia,:,:))';
    eval(['csvwrite(''' name '_CritElOut.out'',ElMax);']);
   
    clear allkk ia kk ElMax
    eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
    end
    eval(['copyfile(''temp.out'',''temp_' name '.out'')']);
    eval(['copyfile(''temp.out'',''' name '.out'')']);
    delete temp.odb temp.out    
    end
 clear name nameout nameodb nameout
 delete temp.odb temp.out abaqus*.*   
end

%%
%---------------------------------------
% POSTPROCESS ALL
%---------------------------------------
left_color = [0. 0. 0.];
right_color = [0. 0. 0.];
%---------------------------------------
% POSTPROCESS UTs
%---------------------------------------
for i=1:nbExps  
name=char(names(i))
eval(['disp('' MESSAGE: ' name ' STARTING POSTPROCESSING'')']);
eval(['copyfile(''temp_' name '.out'',''temp.out'')']);
% eval(['copyfile(''' name '.out'',''temp.out'')']); % IF ALREADY COPIED AND TEMP DELETED
%---------------------------------------
% UTs
%---------------------------------------
if (strncmp(name,'UT_',3)==1)
    eval(['EXP =  experiments.' name ';']); 
    u_beg = 0.002;
    u_end = EXP(length(EXP(:,1)):length(EXP(:,1)),1);
    [kk, ia, ic]  = unique(EXP(:,1),'first');
    ia=sort(ia);
    EXP=EXP(ia,:);
    %
    figure(1); 
    set(1,'defaultAxesColorOrder',[left_color; right_color]);
    plot(EXP(:,1), EXP(:,2), '-g'); hold on;
    eval(['text(EXP(end,1),EXP(end,2),''\leftarrow ' name(4:5) ''',''Color'',''k'')']);
    %---------------------------------------
    [B, t_sim, s] = readColData('temp.out',5,3);
    [kk, ia, ic]  = unique(s(:,1));
    ia=sort(ia);
    s = s(ia,:);
    displ = linspace(u_beg,u_end)';
    for j=1:75, coeff(j)=1;end
    for j=76:95, coeff(j)=1;end
    for j=96:100, coeff(j)=1;end
    %Force pchip
    %F_FEA=pchip(s(:,1),smooth(s(:,2)),displ); %EQ vM
    %plot(displ*6., F_FEA/100., '--b');
    %clear F_FEA; 
    F_FEA=pchip(s(:,1),abs(s(:,4)),displ); %EQ MISES
    plot(displ, F_FEA, '-r'); hold on;
    eval(['text(displ(end),F_FEA(end),''\leftarrow ' name(4:5) ''',''Color'',''r'')']);
    F_exp=pchip(EXP(:,1),abs(EXP(:,2)),displ);
    plot(displ, F_exp, '--k'); hold on;
    
    err(i) = sum(coeff'.*(F_FEA-F_exp).^2)/10000.*0.2; %LESS WEIGHT ON UTs
    eval(['disp('' MESSAGE: ' name ' POSTPROCESSING COMPLETED'')']);
    clear EXP F_FEA F_exp displ s B t_sim ia ic j kk u_beg u_end coeff name; 
%---------------------------------------
% NT20
%---------------------------------------    
elseif (strncmp(name,'NT20_',5)==1)
    eval(['EXP =  experiments.' name ';']); 
    %CHECK UNIQUE U
    [kk, ia, ic]  = unique(EXP(:,2),'first');
    ia=sort(ia);
    EXP=EXP(ia,:);
    %
    u_beg = 0.001;
    u_end = EXP(end,2);
    figure(2);
    set(2,'defaultAxesColorOrder',[left_color; right_color]);
    yyaxis left
    plot(EXP(:,2), EXP(:,3), '-g'); hold on;
    eval(['text(EXP(end,2),EXP(end,3),''\leftarrow NT20-' name(6:7) ''',''Color'',''k'')']); 
    yyaxis right
    plot(EXP(:,2), EXP(:,4)+0., '-g'); hold on;
    eval(['text(EXP(end,2),EXP(end,4),''\leftarrow NT20-' name(6:7) ''',''Color'',''k'')']); 
    %---------------------------------------
    [B, t_sim, s] = readColData('temp.out',22,2);
    [kk, ia, ic]  = unique(s(:,1));
    ia=sort(ia);
    s = s(ia,:);
    displ = linspace(u_beg,u_end)';
	%-----------------------------------
	%Interpolate Experiments
	F_exp=pchip((EXP(:,2)),(EXP(:,3)),displ);
	E_exp=pchip((EXP(:,2)),(EXP(:,4)),displ);
	yyaxis left
	plot(displ, F_exp+0., '--k'); 
	yyaxis right
	plot(displ, E_exp+0., '--k'); hold on;
	%-----------------------------------
    for j=1:75, coeff(j)=1;end
    for j=76:95, coeff(j)=1;end
    for j=96:100, coeff(j)=1;end
    %-----------------------------------
    %SHIFT SIM DISPL TO MATCH ELASTIC RANGE OF EXP
    eval(['displ=displ+deltaU.' name ';']);
    %-----------------------------------	
	%Interpolate SIMULATIONS
    eval(['F_FEA=pchip((2*s(:,18))+deltaU.' name ',(-4./1000.*s(:,19)),displ);']);
    eval(['E_FEA=pchip((2*s(:,18))+deltaU.' name ',log(s(:,21)./(1.0/scalfact)+1.),displ);']);
	%-----------------------------------	
    yyaxis left
    plot(displ, F_FEA+0., '-r'); hold on;
    eval(['text(displ(end),F_FEA(end),''\leftarrow NT20-' name(6:7) ''',''Color'',''r'')']);
    yyaxis right
    plot(displ, E_FEA+0., '-r'); hold on;
    eval(['text(displ(end),E_FEA(end),''\leftarrow NT20-' name(6:7) ''',''Color'',''r'')']);
    % PLOT ALL FORCE
    figure(10);
	plot(displ, F_FEA+0., '-r'); hold on;
	eval(['displ=displ-deltaU.' name ';']);
    plot(displ, F_exp+0., '-k'); hold on;
	eval(['text(EXP(end,2),EXP(end,3),''\leftarrow NT20-' name(6:7) ''',''Color'',''k'')']);
    eval(['text(displ(end),F_FEA(end),''\leftarrow NT20-' name(6:7) ''',''Color'',''r'')']);

    % TRIAX AND LODE
    eval(['Triax.' name '=s(:,13);']); 
    eval(['Lode.' name '=s(:,14);']);
    eval(['PEEQ.' name '=s(:,15);']);
    
    figure(20); hold on;%PLOT TRIAX
    plot(s(:,13), s(:,15), '-k'); 
    eval(['text(s(end,13),s(end,15),''\leftarrow NT20-' name(6:7) ''',''Color'',''k'')']);   
    
    figure(21); hold on;%PLOT LODE
    plot(s(:,14), s(:,15), '-k'); 
    eval(['text(s(end,14),s(end,15),''\leftarrow NT20-' name(6:7) ''',''Color'',''k'')']);   
       
    err(i) = sum(coeff'.*(F_FEA-F_exp).^2)%;+100.*sum(coeff'.*(E_FEA-E_exp).^2);
    eval(['disp('' MESSAGE: ' name ' POSTPROCESSING COMPLETED'')']);
    clear EXP F_FEA F_exp E_FEA E_exp displ s B t_sim ia ic j kk u_beg u_end coeff name; 
%---------------------------------------
% NT6
%---------------------------------------    
elseif (strncmp(name,'NT6_',3)==1)
    eval(['EXP =  experiments.' name ';']); 
    %CHECK UNIQUE U
    [kk, ia, ic]  = unique(EXP(:,2),'first');
    ia=sort(ia);
    EXP=EXP(ia,:);
    %
    u_beg = 0.001;
    u_end = EXP(end,2);
    figure(3);
    set(3,'defaultAxesColorOrder',[left_color; right_color]);
    yyaxis left
    plot(EXP(:,2), EXP(:,3), '-g'); hold on;
    eval(['text(EXP(end,2),EXP(end,3),''\leftarrow NT6-' name(5:6) ''',''Color'',''k'')']); 
    yyaxis right
    plot(EXP(:,2), EXP(:,4)+0., '-g'); hold on;
    eval(['text(EXP(end,2),EXP(end,4),''\leftarrow NT6-' name(5:6) ''',''Color'',''k'')']); 
    %---------------------------------------
    [B, t_sim, s] = readColData('temp.out',22,2);
    [kk, ia, ic]  = unique(s(:,1));
    ia=sort(ia);
    s = s(ia,:);
    displ = linspace(u_beg,u_end)';
 	%-----------------------------------
	%Interpolate Experiments
	F_exp=pchip((EXP(:,2)),(EXP(:,3)),displ);
	E_exp=pchip((EXP(:,2)),(EXP(:,4)),displ);
    yyaxis left
	plot(displ, F_exp+0., '--k'); 
    yyaxis right
	plot(displ, E_exp+0., '--k'); hold on;
	%-----------------------------------
    for j=1:75, coeff(j)=1;end
    for j=76:95, coeff(j)=1;end
    for j=96:100, coeff(j)=1;end
    %-----------------------------------
    %SHIFT SIM DISPL TO MATCH ELASTIC RANGE OF EXP
    eval(['displ=displ+deltaU.' name ';']);
    %-----------------------------------	
	%Interpolate SIMULATIONS
    eval(['F_FEA=pchip((2*s(:,18))+deltaU.' name ',(-4./1000.*s(:,19)),displ);']);
    eval(['E_FEA=pchip((2*s(:,18))+deltaU.' name ',log(s(:,21)./(1.0/scalfact)+1.),displ);']);
	%-----------------------------------	
    yyaxis left
    plot(displ, F_FEA+0., '-r'); hold on;
    eval(['text(displ(end),F_FEA(end),''\leftarrow NT6-' name(5:6) ''',''Color'',''r'')']);
    yyaxis right
    plot(displ, E_FEA+0., '-r'); hold on;
    eval(['text(displ(end),E_FEA(end),''\leftarrow NT6-' name(5:6) ''',''Color'',''r'')']);
    % PLOT ALL FORCE
    figure(10);
	plot(displ, F_FEA+0., '-r'); hold on;
	eval(['displ=displ-deltaU.' name ';']);
    plot(displ, F_exp+0., '-k'); hold on;
	eval(['text(EXP(end,2),EXP(end,3),''\leftarrow NT6-' name(5:6) ''',''Color'',''k'')']);
    eval(['text(displ(end),F_FEA(end),''\leftarrow NT6-' name(5:6) ''',''Color'',''r'')']);

    % TRIAX AND LODE
    eval(['Triax.' name '=s(:,13);']); 
    eval(['Lode.' name '=s(:,14);']);
    eval(['PEEQ.' name '=s(:,15);']);
    
    figure(20); hold on;%PLOT TRIAX
    plot(s(:,13), s(:,15), '-k'); 
    eval(['text(s(end,13),s(end,15),''\leftarrow NT6-' name(5:6) ''',''Color'',''k'')']);   
    
    figure(21); hold on;%PLOT LODE
    plot(s(:,14), s(:,15), '-k'); 
    eval(['text(s(end,14),s(end,15),''\leftarrow NT6-' name(5:6) ''',''Color'',''k'')']);   
       
    err(i) = sum(coeff'.*(F_FEA-F_exp).^2)%;+100.*sum(coeff'.*(E_FEA-E_exp).^2);
    eval(['disp('' MESSAGE: ' name ' POSTPROCESSING COMPLETED'')']);
    clear EXP F_FEA F_exp E_FEA E_exp displ s B t_sim ia ic j kk u_beg u_end coeff name; 
%---------------------------------------
%CENTRAL HOLE CH
%---------------------------------------    
elseif (strncmp(name,'CH_',3)==1)
    eval(['EXP =  experiments.' name ';']); 
    %CHECK UNIQUE U
    [kk, ia, ic]  = unique(EXP(:,2),'first');
    ia=sort(ia);
    EXP=EXP(ia,:);
    %
    u_beg = 0.001;
    u_end = EXP(end,2);
    figure(4);
    set(4,'defaultAxesColorOrder',[left_color; right_color]);
    yyaxis left
    plot(EXP(:,2), EXP(:,3), '-g'); hold on;
    eval(['text(EXP(end,2),EXP(end,3),''\leftarrow CH-' name(4:5) ''',''Color'',''k'')']); 
    yyaxis right
    plot(EXP(:,2), EXP(:,4)+0., '-g'); hold on;
    plot(EXP(:,2), EXP(:,5)+0., '-g'); hold on;
    eval(['text(EXP(end,2),EXP(end,4),''\leftarrow CH-' name(4:5) ''',''Color'',''k'')']); 
    eval(['text(EXP(end,2),EXP(end,5),''\leftarrow CH-' name(4:5) ''',''Color'',''k'')']);
    %---------------------------------------
    [B, t_sim, s] = readColData('temp.out',22,2);
    [kk, ia, ic]  = unique(s(:,1));
    ia=sort(ia);
    s = s(ia,:);
    displ = linspace(u_beg,u_end)';
		%-----------------------------------
	%Interpolate Experiments
	F_exp=pchip((EXP(:,2)),(EXP(:,3)),displ);
	E_exp=pchip((EXP(:,2)),(EXP(:,4)),displ);
	yyaxis left
	plot(displ, F_exp+0., '--k'); 
	yyaxis right
	plot(displ, E_exp+0., '--k'); hold on;
	%-----------------------------------
    for j=1:75, coeff(j)=1;end
    for j=76:95, coeff(j)=1;end
    for j=96:100, coeff(j)=1;end
    %-----------------------------------
    %SHIFT SIM DISPL TO MATCH ELASTIC RANGE OF EXP
    eval(['displ=displ+deltaU.' name ';']);
    %-----------------------------------	
	%Interpolate SIMULATIONS
    eval(['F_FEA=pchip((2*s(:,18))+deltaU.' name ',(-4./1000.*s(:,19)),displ);']);
    eval(['E_FEA=pchip((2*s(:,18))+deltaU.' name ',log(s(:,21)./(1.0/scalfact)+1.),displ);']);
	%-----------------------------------	
    yyaxis left
    plot(displ, F_FEA+0., '-r'); hold on;
    eval(['text(displ(end),F_FEA(end),''\leftarrow CH-' name(4:5) ''',''Color'',''r'')']);
    yyaxis right
    plot(displ, E_FEA+0., '-r'); hold on;
    eval(['text(displ(end),E_FEA(end),''\leftarrow CH-' name(4:5) ''',''Color'',''r'')']);
    % PLOT ALL FORCE
    figure(10);
	plot(displ, F_FEA+0., '-r'); hold on;
	eval(['displ=displ-deltaU.' name ';']);
    plot(displ, F_exp+0., '-k'); hold on;
	eval(['text(EXP(end,2),EXP(end,3),''\leftarrow CH-' name(4:5) ''',''Color'',''k'')']);
    eval(['text(displ(end),F_FEA(end),''\leftarrow CH-' name(4:5) ''',''Color'',''r'')']);

    % TRIAX AND LODE
    eval(['Triax.' name '=s(:,13);']); 
    eval(['Lode.' name '=s(:,14);']);
    eval(['PEEQ.' name '=s(:,15);']);
    
    figure(20); hold on;%PLOT TRIAX
    plot(s(:,13), s(:,15), '-k'); 
    eval(['text(s(end,13),s(end,15),''\leftarrow CH-' name(4:5) ''',''Color'',''k'')']);   
    
    figure(21); hold on;%PLOT LODE
    plot(s(:,14), s(:,15), '-k'); 
    eval(['text(s(end,14),s(end,15),''\leftarrow CH-' name(4:5) ''',''Color'',''k'')']);   
  
    err(i) = sum(coeff'.*(F_FEA-F_exp).^2)%;+100.*sum(coeff'.*(E_FEA-E_exp).^2);
    eval(['disp('' MESSAGE: ' name ' POSTPROCESSING COMPLETED'')']);
    clear EXP F_FEA F_exp E_FEA E_exp displ s B t_sim ia ic j kk u_beg u_end coeff name; 
%---------------------------------------
% SHEAR
%---------------------------------------    
elseif (strncmp(name,'SH_',3)==1)
    eval(['EXP =  experiments.' name ';']); 
    %CHECK UNIQUE U
    [kk, ia, ic]  = unique(EXP(:,2),'first');
    ia=sort(ia);
    EXP=EXP(ia,:);
    %
    u_beg = 0.001;
    u_end = EXP(end,2);
    figure(5);
    set(5,'defaultAxesColorOrder',[left_color; right_color]);
    plot(EXP(:,2), EXP(:,3), '-g'); hold on;
    eval(['text(EXP(end,2),EXP(end,3),''\leftarrow SH-' name(4:5) ''',''Color'',''k'')']);  
    %---------------------------------------
    [B, t_sim, s] = readColData('temp.out',3,3);
    [kk, ia, ic]  = unique(s(:,1));
    ia=sort(ia);
    s = s(ia,:);
    displ = linspace(u_beg,u_end)';
	 	%-----------------------------------
	%Interpolate Experiments
	F_exp=pchip((EXP(:,2)),(EXP(:,3)),displ);
	plot(displ, F_exp+0., '--k'); 
	%-----------------------------------
    for j=1:75, coeff(j)=1;end
    for j=76:95, coeff(j)=1;end
    for j=96:100, coeff(j)=1;end
    %-----------------------------------
    %SHIFT SIM DISPL TO MATCH ELASTIC RANGE OF EXP
    eval(['displ=displ+deltaU.' name ';']);
    %-----------------------------------	
	%Interpolate SIMULATIONS
    F_FEA=pchip((s(:,1)),(-2./1000.*s(:,2)),displ);
	%-----------------------------------	
    plot(displ, F_FEA+0., '-r'); hold on;
    eval(['text(displ(end),F_FEA(end),''\leftarrow SH-' name(4:5) ''',''Color'',''r'')']);
    % PLOT ALL FORCE
    figure(10);
	plot(displ, F_FEA+0., '-r'); hold on;
	eval(['displ=displ-deltaU.' name ';']);
    plot(displ, F_exp+0., '-k'); hold on;
	eval(['text(EXP(end,2),EXP(end,3),''\leftarrow SH-' name(4:5) ''',''Color'',''k'')']);
    eval(['text(displ(end),F_FEA(end),''\leftarrow SH-' name(4:5) ''',''Color'',''r'')']);
    
    eval(['ElMax=csvread(''' name '_CritElOut.out'');']);
    
    eval(['Triax.' name '=ElMax(:,3);']); 
    eval(['Lode.' name '=ElMax(:,4);']); 
    eval(['PEEQ.' name '=ElMax(:,5);']); 
    
    figure(20); hold on;%PLOT TRIAX
    eval(['plot(Triax.' name ',PEEQ.' name ',''-k'');']);
    eval(['text(Triax.' name '(end),PEEQ.' name '(end),''\leftarrow SH-' name(4:5) ''',''Color'',''k'')']);   
    
    figure(21); hold on;%PLOT LODE
    eval(['plot(Lode.' name ',PEEQ.' name ',''-k'');']);
    eval(['text(Lode.' name '(end),PEEQ.' name '(end),''\leftarrow SH-' name(4:5) ''',''Color'',''k'')']);   
        
    err(i) = 4.*sum(coeff'.*(F_FEA-F_exp).^2);%+100.*sum(coeff'.*(E_FEA-E_exp).^2);
    eval(['disp('' MESSAGE: ' name ' POSTPROCESSING COMPLETED'')']);
    clear EXP F_FEA F_exp E_FEA E_exp displ s B t_sim ia ic j kk u_beg u_end coeff name NODES ELEMENTS centerEls ElMax;       
end
delete temp.out temp.odb
end
%--------------------------------------------------------------
% SAVE TRIAX< LODE PEEQ TO FILE / GET MAX
%--------------------------------------------------------------
pathnames=(fieldnames(PEEQ))
ll=length(pathnames);
for j=1:ll
    clear aaa
    eval(['aaa(:,1)=Triax.' char(pathnames(j)) ';']);
    eval(['aaa(:,2)=Lode.' char(pathnames(j)) ';']);   
    eval(['aaa(:,3)=PEEQ.' char(pathnames(j)) ';']);
    eval(['csvwrite(''_path_' char(pathnames(j)) '.csv'',aaa);']);
     
    eval(['ymax(j)=max(PEEQ.' char(pathnames(j)) ');']);
    clear aaa
end
ymax=1.05*max(ymax);
%--------------------------------------------------------------
% SAVE PLOTS
%--------------------------------------------------------------
if ishandle(1)
figure(1);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
xlabel('Plastic Strain [-]','FontSize',12);
ylabel('True Stress [MPa]','FontSize',12);yaxismin = 10000;
yaxismax = 0;
for i=1:nbExps
    name=char(names(i));
    if (strncmp(name,'UT_',3)==1)
    eval(['EXP =  experiments.' name ';']); 
    yaxismin = min(yaxismin,min(EXP(:,2)));
    yaxismax = max(yaxismax,max(EXP(:,2)));
    clear EXP name
    end
end
axis([0 inf yaxismin*.95 yaxismax*1.05])
axis square;
box on
saveas(gcf,['_final_UT.png'],'png');
savefig(1,'_final_UT');
hold off;
clear yaxismin yaxismax
end
 %--------------------------------------------------------------
if ishandle(2)
figure(2);
yyaxis left
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
xlabel('Displacement [mm]','FontSize',12);
ylabel('Force [kN]','FontSize',12);
xaxismax = 0;
yaxis1max = 0;
yaxis2max = 0;
for i=1:nbExps
    name=char(names(i));
    if (strncmp(name,'NT20_',5)==1)
    eval(['EXP =  experiments.' name ';']); 
    xaxismax = max(xaxismax,max(EXP(:,2)));
    yaxis1max = max(yaxis1max,max(EXP(:,3)));
    yaxis2max = max(yaxis2max,max(EXP(:,4)));
    clear EXP name
    end
end
axis([0 xaxismax*1.02 0 yaxis1max*1.05])
box on
yyaxis right
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
% xlabel('Displacement [mm]','FontSize',12);
ylabel('Log. Ax. Strain [-]','FontSize',12);
axis([0 xaxismax*1.02 0 yaxis2max*1.5])
box on
axis square
saveas(gcf,['_final_NT20.png'],'png');
savefig(2,'_final_NT20');
hold off;
clear xaxismax yaxis1max yaxis2max
end
%--------------------------------------------------------------
if ishandle(3)
figure(3);
yyaxis left
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
xlabel('Displacement [mm]','FontSize',12);
ylabel('Force [kN]','FontSize',12);
xaxismax = 0;
yaxis1max = 0;
yaxis2max = 0;
for i=1:nbExps
    name=char(names(i));
    if (strncmp(name,'NT6_',4)==1)
    eval(['EXP =  experiments.' name ';']); 
    xaxismax = max(xaxismax,max(EXP(:,2)));
    yaxis1max = max(yaxis1max,max(EXP(:,3)));
    yaxis2max = max(yaxis2max,max(EXP(:,4)));
    clear EXP name
    end
end
axis([0 xaxismax*1.02 0 yaxis1max*1.05])
box on
yyaxis right
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
% xlabel('Displacement [mm]','FontSize',12);
ylabel('Log. Ax. Strain [-]','FontSize',12);
axis([0 xaxismax*1.02 0 yaxis2max*1.5])
axis square;
box on
saveas(gcf,['_final_NT6.png'],'png');
savefig(3,'_final_NT6');
hold off;
clear xaxismax yaxis1max yaxis2max
end
%--------------------------------------------------------------
if ishandle(4)
figure(4);
yyaxis left
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
xlabel('Displacement [mm]','FontSize',12);
ylabel('Force [kN]','FontSize',12);
xaxismax = 0;
yaxis1max = 0;
yaxis2max = 0;
for i=1:nbExps
    name=char(names(i));
    if (strncmp(name,'CH_',3)==1)
    eval(['EXP =  experiments.' name ';']); 
    xaxismax = max(xaxismax,max(EXP(:,2)));
    yaxis1max = max(yaxis1max,max(EXP(:,3)));
    yaxis2max = max(yaxis2max,max(max(EXP(:,5),max(EXP(:,4)))));
    clear EXP name
    end
end
axis([0 xaxismax*1.02 0 yaxis1max*1.05])
box on
yyaxis right
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
% xlabel('Displacement [mm]','FontSize',12);
ylabel('Log. Ax. Strain [-]','FontSize',12);
axis([0 xaxismax*1.02 0 yaxis2max*1.5])
axis square;
box on
saveas(gcf,['_final_CH.png'],'png');
savefig(4,'_final_CH');
hold off;
clear xaxismax yaxis1max yaxis2max
end
%--------------------------------------------------------------
if ishandle(5)
figure(5);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
xlabel('Displacement [mm]','FontSize',12);
ylabel('Force [kN]','FontSize',12);
xaxismax = 0;
yaxis1max = 0;
yaxis2max = 0;
for i=1:nbExps
    name=char(names(i));
    if (strncmp(name,'SH_',3)==1)
    eval(['EXP =  experiments.' name ';']); 
    xaxismax = max(xaxismax,max(EXP(:,2)));
    yaxis1max = max(yaxis1max,max(EXP(:,3)));
    clear EXP name
    end
end
axis([0 xaxismax*1.02 0 yaxis1max*1.05])
axis square;
box on
saveas(gcf,['_final_SH.png'],'png');
savefig(5,'_final_SH');
hold off;
end
%--------------------------------------------------------------
if ishandle(10)
figure(10);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
xlabel('Displacement [mm]','FontSize',12);
ylabel('Force [kN]','FontSize',12);
axis([0 inf 0 inf])
axis square;
box on
saveas(gcf,['_final_AllFrac.png'],'png');
savefig(10,'_final_AllFrac');
hold off;
end
%--------------------------------------------------------------
if ishandle(20)
figure(20);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
xlabel('Triaxiality [-]','FontSize',12);
ylabel('Eq. Plastic Strain [-]','FontSize',12);
axis([0.0 inf 0.0 ymax])
set(gca,'XTick',[-0.33,0,0.1,0.2,0.33,0.45,0.57,0.67,0.8,0.9,1]);
axis square;
box on
saveas(gcf,['_final_Triax_PEEQ.png'],'png');
savefig(21,'_final_Triax_PEEQ');
hold off;
end
%--------------------------------------------------------------
if ishandle(21)
figure(21);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
xlabel('Lode Angle Parameter [-]','FontSize',12);
ylabel('Eq. Plastic Strain [-]','FontSize',12);
axis([-1.0 1.0 0.0 ymax])
axis square;
box on
saveas(gcf,['_final_Lode_PEEQ.png'],'png');
savefig(21,'_final_Lode_PEEQ');
hold off;
end
%------------------------------------------------
% FINISH POSTPROCESSING CLEAN UP
%------------------------------------------------
delete temp*.*
%------------------------------------------------
% Calculate EMC
%------------------------------------------------
EMC_main_ident
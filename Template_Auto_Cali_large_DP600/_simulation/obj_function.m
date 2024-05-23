function result = obj_function(pv)
global  experiments nbExps deltaU YLD2000_pre Hill48_pre minres Loop;
global a0 P Pi Sigma K0;
delete screen.txt;
diary screen.txt;
scalfact = 1; % Scaling compared to the bigest set of the 3 specimens (big: 1, medium: 2, small:4)
scalfact = scalfact*5;
%
%*DENSITY
%7.8e-6
%*USER MATERIAL,CONSTANTS=46
%**Physical constants-------------------------------------------------------------------------------
%**       E          v         Cp         Xi          -          -          -  Fflag(0/1) (BE/FE)   
%     210.0,       0.3,     500.0,       0.0,       0.0,       0.0,       0.0,       0.0, 
%**Yield function and flow rule>HILL-NAFR-----------------------------------------------------------
%**     P12,       P22,       P33,       G12,       G22,       G33,
%      -0.5,       1.0,       3.0,      -0.5,       1.0,       3.0,       0.0,       0.0,
%**Hardening----------------------------------------------------------------------------------------
%**       A         B          n         gs           -          -          -       (1) HARD JCX    
%**      s0        Q1         C1         Q2          C2          -          -       (2) HARD VOCE   
%**       A        e0          n         s0          Q1          C1       alpha        (3) HARD MSV    
%     0.630,     0.380,     9.000,     8.300,     0.040,     0.000,       0.0,        3,
%**SRH+TS>JC----------------------------------------------------------------------------------------
%**  eps0dot         C          m         T0         Tr         Tm    epsAdot         -
%    5.0e-7,     0.000,    0.8350,      25.0,      25.0,    1500.0,       0.0,      0.0,
%**Failure criterion -> ADDFail (+) Tc, (-) EQPSmax-------------------------------------------------
%**     Wcr          -         -          -           -          -    ADDFail      (0/1) FAILURE CL 
%**       a          b         c          n          D4         D5    ADDFail      (2/3) FAILURE HC 
%**      D1         D2        D3         gf          D4         D5    ADDFail      (4/5) FAILURE JCX
%       1.0,      0.0,       0.0,        0.0,       0.0,       0.0,       0.0,        2,
%**Post initiation damage---------------------------------------------------------------------------
%**      D0,        Dc,        mD,     DcMax,      bMinS,       k,
%       1.0,       2.0,      0.00,      1.00,        0.1,     1.0,
%Evaluation up to
All_x_Fmax ='All'; %'All' or 'Fmax' 

% READ DATA FROM PREP FILE
Angle    = Hill48_pre(1,:)
SigY     = Hill48_pre(2,:)
SigRatio = Hill48_pre(3,:)
RRatio   = Hill48_pre(4,:)

Swift_A  = Hill48_pre(7,1); %
Swift_e0 = Hill48_pre(7,2);%
Swift_n  = Hill48_pre(7,3);%

Voce_S0  = Hill48_pre(7,4);%
Voce_Q1  = Hill48_pre(7,5);%
Voce_C1  = Hill48_pre(7,6);%

SValpha =  Hill48_pre(7,7);%

YS0  = SigY(1);
YS45 = SigY(4);
YS90 = SigY(7);
YSB  = SigY(1);
% YSB  = pv(1).*SigY(1);

% CALC Ps and Gs

p12 = Hill48_pre(5,1);%0.5*(pv(1)).^2-((1./SigRatio(7))^2-1);%-pv(1) , 
p22 = Hill48_pre(5,2);%(1./SigRatio(7)).^2;% pv(2)  
p44 = Hill48_pre(5,3);%(2.*(1./SigRatio(4))).^2-pv(1).^2;%pv(3) 
% 
% p12 = 0.5*((YS0./YSB).^2-(YS0./YS90).^2-1);%-pv(1) , 
% p22 = (YS0./YS90).^2;% pv(2)  
% p44 = (2.*YS0./YS45).^2-(YS0./YSB).^2;%pv(3)  

% p12 = 0.5*((YS0./YSB).^2-(YS0./YS90).^2-1.);% -pv(1);%Hill48_pre(5,1);%
% p22 = (YS0./YS90).^2;%  pv(2);%Hill48_pre(5,2);%
% p44 = (2.*YS0./YS45).^2-(YS0./YSB).^2;% pv(3);%Hill48_pre(5,3);%

g12 = Hill48_pre(5,4) ;%-r_00./(1+r_00); %p12;%-pv(1); 
g22 = Hill48_pre(5,5) ;%r_00./r_90.*(1+r_90)./(1+r_00); % p22;%  pv(2);%
g44 = Hill48_pre(5,6) ;%(1.+2.*r_45)./r_90.*(r_00+r_90)./(1+r_00); %p44;% pv(3);%

spectemp = [25.0];
Etemp    = [200.E3];
for iii=1:length(spectemp)
eval(['name = ''material.dat'';']);
fid=fopen(name,'w');
fprintf(fid,'*User Material, constants=46 \n');
fprintf(fid,['' num2str(Etemp(iii)) ', 0.3, 4.49E8, 0.9, 0.0, 0.0, 0.0, 0.0,\n']); 
fprintf(fid,['' mat2str(p12) ', ' mat2str(p22) ', ' mat2str(p44) ', ' mat2str(g12) ', ' mat2str(g22) ', ' mat2str(g44) ', 0.0, 0.0,\n']);
fprintf(fid,['' mat2str(Swift_A) ', ' mat2str(Swift_e0) ', ' mat2str(Swift_n) ', ' mat2str(Voce_S0) ', ' mat2str(Voce_Q1) ', ' mat2str(Voce_C1) ', ' mat2str(SValpha) ', 3.0,\n']);
fprintf(fid,['1.E-3, 0.0, 0.0 , ' num2str(spectemp(iii)) ', 25.0, 1200., 20., 0.0,\n']);
% fprintf(fid,['0.001, 0.0, ' mat2str(pv(8)) ', ' num2str(spectemp(iii)) ', ' mat2str(pv(9)) 'E2, ' mat2str(pv(10)) 'E2, 2.0, 0.0,\n']);
% fprintf(fid,'10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\n');
fprintf(fid,'2.0, 10.0, 0.1, 0.1, 0.0, 0.0, 0.0, 2.0,\n');
fprintf(fid,'1.0, 2.0, 0.00, 1.0, 0.1, 1.0,\n');
fprintf(fid,'*DENSITY \n');
fprintf(fid,'7.85E-9,');    
fclose(fid);
end
% ---------------------------------
% CREATE Y and R plot
% ---------------------------------
adeg = [0:1:90];
sig_r=1./(sqrt(cos(adeg.*pi()./180).^4+sin(adeg.*pi()./180).^2.*cos(adeg.*pi()./180).^2.*(2.*p12+p44)+sin(adeg.*pi()./180).^4.*p22));
r=(((g44+2.*g12-g22-1).*sin(adeg.*pi()./180).^2.*cos(adeg.*pi()./180).^2-g12)./((1-g22).*cos(adeg.*pi()./180).^2+g12+g22));
fig=openfig('__Y-R_ratios_Angle.fig');

allaxes = findall(fig, 'type', 'axes');
if (strncmp(get(allaxes(1),'YAxisLocation'),'left',4)==1) 
    set(fig,'CurrentAxes',allaxes(1));
    plot(adeg,sig_r,'-k','linewidth',1.0);
    set(fig,'CurrentAxes',allaxes(2));
    plot(adeg,r,'-b','linewidth',1.0);
elseif (strncmp(get(allaxes(1),'YAxisLocation'),'right',5)==1)
    set(fig,'CurrentAxes',allaxes(2));
    plot(adeg,sig_r,'-k','linewidth',1.0);
    set(fig,'CurrentAxes',allaxes(1));
    plot(adeg,r,'-b','linewidth',1.0);
end

annotation('textbox',...
    [0.2 0.85 0.15 0.05],...
    'String',{'Hill''48'},...
    'FontSize',18,...
    'FontName','Tahoma',...
    'Color',[0 0 0]);

% % axis square
saveas(gcf,['_curve_' mat2str(Loop) '_Y-R_ratios_Angle.png'],'png');
% 
eval(['movefile(''_curve_' mat2str(Loop) '_Y-R_ratios_Angle.png'',''_overview/_iterations'')']);
% 
clear sig_r r adeg 
close all
%---------------------------------------
% Cluster folder cleanup commands
%---------------------------------------
delete lsf*.*;
delete abaqus*.*;
delete temp*.*;
%---------------------------------------
% COPY SIMs to TEMPSIM files % LAUNCH COMPUTATIONS
%---------------------------------------
launchcount = 0;
names=fieldnames(experiments);
for i=1:nbExps
name=char(names(i))
%---------------------------------------
% UTs
%---------------------------------------
if (strncmp(name,'UT_',3)==1)
  eval(['copyfile(''' name '.inp'',''temp_' name '.inp'')']);
  %exe_command = ['!bsub -n 1 -W 2:00 -R "ib rusage[mem=500,scratch=500]" "abaqus job=temp_' name ' user=VUMAT_explicit.f double cpus=1 memory="500mb" scratch=\$TMPDIR"'];
  exe_command = ['!sbatch -n 1 -t 0-2 --mem-per-cpu=500 --tmp=500 --wrap "abaqus job=temp_' name ' user=VUMAT_explicit.f double cpus=1 scratch=\$TMPDIR"'];
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
  WIDTH =  tempdata(1,8);
  eval(['DISPL=tempdata(end,2)-deltaU.' name '+0.0001;']);
  DTIME=tempdata(end,1); 
  TRAMP1= 0.01*DTIME
  RAMP1= 0.001*DISPL  
  eval(['fid=fopen(''' name '_paras.inp'',''wt'');']);
  %
  DT	 = DTIME/1.E5;
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
  %exe_command = ['!bsub -n 12 -W 2:00 -R "ib rusage[mem=300,scratch=300]" "abaqus job=temp_' name ' user=VUMAT_explicit.f double cpus=12 memory="3600mb" scratch=\$TMPDIR"'];
  exe_command = ['!sbatch -n 12 0-2 --mem-per-cpu=300 --tmp=300 --wrap "abaqus job=temp_' name ' user=VUMAT_explicit.f double cpus=12 memory="3600mb" scratch=\$TMPDIR"'];
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
  
  eval(['DISPL=tempdata(end,2)-deltaU.' name '+0.0001;']);
  DTIME=tempdata(end,1);
  TRAMP1= 0.01*DTIME
  RAMP1= 0.001*DISPL  
  eval(['fid=fopen(''' name '_paras.inp'',''wt'');']);
  %
  if abs(tempdata(1,8))<=6.5
   copyfile('CH_r2_5.geo','CH.geo','f');
   copyfile('CH_r2_5.py','CH.py','f');
 %  WIDTH = 5./5.;%tempdata(1,8)/5.;
   WIDTH = tempdata(1,7)/20.;
  end
  if abs(tempdata(1,8))>6.5
   copyfile('CH_r4.geo','CH.geo','f');
   copyfile('CH_r4.py','CH.py','f');
   WIDTH = 8./8.;%tempdata(1,8)/8.;
  end
  
  DT	 = DTIME/1.E5;
  fprintf(fid,'*PARAMETER\n');
  fprintf(fid,'DISPL = %f\n',DISPL/2.);
  fprintf(fid,'THICKNESS = %f\n',THICKNESS);
  fprintf(fid,'WIDTH = %f\n',WIDTH);
  fprintf(fid,'DTIME = %f\n',DTIME);
  fprintf(fid,'DT = %f\n',DT);
  fprintf(fid,'TRAMP1 = %f\n',TRAMP1);
  fprintf(fid,'RAMP1 = %f\n',RAMP1);
  fclose(fid); 
  %
  exe_command = ['!sbatch -n 12 0-2 --mem-per-cpu=300 --tmp=300 --wrap "abaqus job=temp_' name ' user=VUMAT_explicit.f double cpus=12 memory="3600mb" scratch=\$TMPDIR"'];
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
  %WIDTH = 20.;
  WIDTH = tempdata(1,7);
  eval(['DISPL=tempdata(end,2)-deltaU.' name '+0.001;']);
  DTIME=tempdata(end,1);
  TRAMP1= 0.01*DTIME
  RAMP1= 0.001*DISPL  
  eval(['fid=fopen(''' name '_paras.inp'',''wt'');']);
  %
  DT	 = DTIME/2.0E5;
  fprintf(fid,'*PARAMETER\n');
  fprintf(fid,'DISPL = %f\n',DISPL);
  fprintf(fid,'THICKNESS = %f\n',THICKNESS/2.*1.25);
  fprintf(fid,'WIDTH = %f\n',WIDTH/20.);
  fprintf(fid,'DTIME = %f\n',DTIME);
  fprintf(fid,'TRAMP1 = %f\n',TRAMP1);
  fprintf(fid,'RAMP1 = %f\n',RAMP1);
  fprintf(fid,'DT = %f\n',DT);
  fclose(fid); 
  %
  exe_command = ['!sbatch -n 12 0-2 --mem-per-cpu=300 --tmp=300 --wrap "abaqus job=temp_' name ' user=VUMAT_explicit.f double cpus=12 memory="3600mb" scratch=\$TMPDIR"'];
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
     pause(5);
     %---------------------------------------
     if (strncmp(name,'UT_',3)==1)  
     !abaqus cae noGUI=UT.py
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'NT20_',5)==1)  
     !abaqus cae noGUI=NT20.py  
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'NT6_',4)==1)  
     !abaqus cae noGUI=NT6.py  
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'CH_',3)==1)  
     !abaqus cae noGUI=CH.py 
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'SH_',3)==1)  
     !abaqus cae noGUI=SH.py 
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
%     eval(['!copy temp.out' ' ' name '.out ']);
     end
     eval(['copyfile(''temp.out'',''temp_' name '.out'')']);
     delete temp.odb temp.out    
     %---------------------------------------
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
     pause(5);
     %---------------------------------------
     if (strncmp(name,'UT_',3)==1)  
     !abaqus cae noGUI=UT.py
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'NT20_',5)==1)  
     !abaqus cae noGUI=NT20.py  
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'NT6_',4)==1)  
     !abaqus cae noGUI=NT6.py  
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'CH_',3)==1)  
     !abaqus cae noGUI=CH.py 
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     elseif (strncmp(name,'SH_',3)==1)  
     !abaqus cae noGUI=SH.py 
     eval(['disp('' MESSAGE: ' name ' PYTHON SCRIPT DONE'')']);
     end
     eval(['copyfile(''temp.out'',''temp_' name '.out'')']);
     delete temp.odb temp.out    
    end
 clear name nameout nameodb nameout
end
pause(5)
%---------------------------------------
% POSTPROCESS UTs
%---------------------------------------
%%
left_color = [0. 0. 0.];
right_color = [0. 0. 0.];
for i=1:nbExps  
name=char(names(i))
eval(['disp('' MESSAGE: ' name ' STARTING POSTPROCESSING'')']);
eval(['copyfile(''temp_' name '.out'',''temp.out'')']);
%---------------------------------------
% UTs
%---------------------------------------
if (strncmp(name,'UT_',3)==1)
    eval(['EXP =  experiments.' name ';']); 
    u_beg = 0.000;
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
    F_FEA=pchip(s(:,1),(s(:,4)),displ); %EQ vM
    plot(displ, F_FEA, '-r'); hold on;
    eval(['text(displ(end),F_FEA(end),''\leftarrow ' name(4:5) ''',''Color'',''r'')']);
    F_exp=pchip((EXP(:,1)),(EXP(:,2)),displ);
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
    %CHOOSE DURATION
    u_beg = 0.00;
    switch All_x_Fmax
    case 'All' 
        u_end = EXP(end,2);
    case 'Fmax'     
        [ic line]=max(EXP(:,3));
        u_end = max(EXP(line,2))  
    end
    figure(2);
    set(2,'defaultAxesColorOrder',[left_color; right_color]);
    yyaxis left
    plot(EXP(:,2), EXP(:,3), '-g'); hold on;
    eval(['text(EXP(end,2),EXP(end,3),''\leftarrow NT20-' name(6:7) ''',''Color'',''k'')']); 
    yyaxis right
    plot(EXP(:,2), EXP(:,4)+0., '-g'); hold on;
    eval(['text(EXP(end,2),EXP(end,4),''\leftarrow NT20-' name(6:7) ''',''Color'',''k'')']); 
    %---------------------------------------
    [B, t_sim, s] = readColData('temp.out',4,3);
    [kk, ia, ic]  = unique(s(:,1));
    ia=sort(ia);
    s = s(ia,:);
    displ = linspace(u_beg,u_end)';
    for j=1:75, coeff(j)=1;end
    for j=76:95, coeff(j)=1;end
    for j=96:100, coeff(j)=1;end
    %-----------------------------------
    % SHIFT SIM DISPL TO MATCH ELASTIC RANGE OF EXP
    eval(['displ=displ+deltaU.' name ';']);
    %-----------------------------------
    %Force pchip
	eval(['F_FEA=pchip((2*s(:,1)+deltaU.' name '),(-4./1000.*s(:,2)),displ)']);
	yyaxis left
    plot(displ, F_FEA+0., '-r'); hold on;
    eval(['text(displ(end),F_FEA(end),''\leftarrow NT20-' name(6:7) ''',''Color'',''r'')']);
    F_exp=pchip(EXP(:,2),abs(EXP(:,3)),displ);
    plot(displ, F_exp+0., '--k'); 
    %Strain pchip
     E_FEA=pchip((2*s(:,1)),log(s(:,3)./(0.50/scalfact)+1.),displ);
    yyaxis right
    plot(displ, E_FEA+0., '-r'); hold on;
    eval(['text(displ(end),E_FEA(end),''\leftarrow NT20-' name(6:7) ''',''Color'',''r'')']);
    try
    E_exp=pchip((EXP(:,2)),(abs(EXP(:,4))),displ);
    plot(displ, E_exp+0., '--k'); hold on;
    end
    % PLOT ALL FORCE
    figure(10);
    plot(displ, F_exp+0., '-k'); hold on;
    plot(displ, F_FEA+0., '-r'); hold on;
    eval(['text(EXP(end,2),EXP(end,3),''\leftarrow NT20-' name(6:7) ''',''Color'',''k'')']);
    eval(['text(displ(end),F_FEA(end),''\leftarrow NT20-' name(6:7) ''',''Color'',''r'')']);
        
    err(i) = sum(coeff'.*(F_FEA-F_exp).^2);%;+100.*sum(coeff'.*(E_FEA-E_exp).^2);
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
    %CHOOSE DURATION
    u_beg = 0.00;
    switch All_x_Fmax
    case 'All' 
        u_end = EXP(end,2);
    case 'Fmax'     
        [ic line]=max(EXP(:,3));
        u_end = max(EXP(line,2))  
    end
    figure(3);
    set(3,'defaultAxesColorOrder',[left_color; right_color]);
    yyaxis left
    plot(EXP(:,2), EXP(:,3), '-g'); hold on;
    eval(['text(EXP(end,2),EXP(end,3),''\leftarrow NT6-' name(5:6) ''',''Color'',''k'')']); 
    yyaxis right
    plot(EXP(:,2), EXP(:,4)+0., '-g'); hold on;
    eval(['text(EXP(end,2),EXP(end,4),''\leftarrow NT6-' name(5:6) ''',''Color'',''k'')']); 
      %---------------------------------------
    [B, t_sim, s] = readColData('temp.out',4,3);
    [kk, ia, ic]  = unique(s(:,1));
    ia=sort(ia);
    s = s(ia,:);
    displ = linspace(u_beg,u_end)';
    for j=1:75, coeff(j)=1;end
    for j=76:95, coeff(j)=1;end
    for j=96:100, coeff(j)=1;end
    %-----------------------------------
    % SHIFT SIM DISPL TO MATCH ELASTIC RANGE OF EXP
    eval(['displ=displ+deltaU.' name ';']);
    %-----------------------------------
    %Force pchip
	eval(['F_FEA=pchip((2*s(:,1)+deltaU.' name '),(-4./1000.*s(:,2)),displ)']);
    yyaxis left
    plot(displ, F_FEA+0., '-r'); hold on;
    eval(['text(displ(end),F_FEA(end),''\leftarrow NT6-' name(5:6) ''',''Color'',''r'')']);
    F_exp=pchip((EXP(:,2)),(EXP(:,3)),displ);
    plot(displ, F_exp+0., '--k'); 
    %Strain pchip
    E_FEA=pchip((2*s(:,1)),log(s(:,3)./(0.50/scalfact)+1.),displ);
    yyaxis right
    plot(displ, E_FEA+0., '-r'); 
    eval(['text(displ(end),E_FEA(end),''\leftarrow NT6-' name(5:6) ''',''Color'',''r'')']);
    try
    E_exp=pchip((EXP(:,2)),(abs(EXP(:,4))),displ);
    plot(displ, E_exp+0., '--k'); hold on;
    end
    % PLOT ALL FORCE
    figure(10);
    plot(displ, F_exp+0., '-k'); hold on;
    plot(displ, F_FEA+0., '-r'); hold on;
    eval(['text(EXP(end,2),EXP(end,3),''\leftarrow NT6-' name(5:6) ''',''Color'',''k'')']); 
    eval(['text(displ(end),F_FEA(end),''\leftarrow NT6-' name(5:6) ''',''Color'',''r'')']);
    
    err(i) = sum(coeff'.*(F_FEA-F_exp).^2);%;+100.*sum(coeff'.*(E_FEA-E_exp).^2);
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
    %CHOOSE DURATION
    u_beg = 0.00;
    switch All_x_Fmax
    case 'All' 
        u_end = EXP(end,2);
    case 'Fmax'     
        [ic line]=max(EXP(:,3));
        u_end = max(EXP(line,2))  
    end
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
    [B, t_sim, s] = readColData('temp.out',4,3);
    [kk, ia, ic]  = unique(s(:,1));
    ia=sort(ia);
    s = s(ia,:);
    displ = linspace(u_beg,u_end)';
    for j=1:75, coeff(j)=1;end
    for j=76:95, coeff(j)=1;end
    for j=96:100, coeff(j)=1;end
    %-----------------------------------
    % SHIFT SIM DISPL TO MATCH ELASTIC RANGE OF EXP
    eval(['displ=displ+deltaU.' name ';']);
    %-----------------------------------
    %Force pchip
	eval(['F_FEA=pchip((2*s(:,1)+deltaU.' name '),(-4./1000.*s(:,2)),displ)']);
    yyaxis left
    plot(displ, F_FEA+0., '-r'); 
    eval(['text(displ(end),F_FEA(end),''\leftarrow CH-' name(4:5) ''',''Color'',''r'')']);
    F_exp=pchip((EXP(:,2)),(EXP(:,3)),displ);
    plot(displ, F_exp+0., '--k'); 
    %Strain pchip
    E_FEA=pchip((2*s(:,1)),log(s(:,3)./(0.50/scalfact)+1.),displ);
    yyaxis right
    plot(displ, E_FEA+0., '-r'); 
    eval(['text(displ(end),E_FEA(end),''\leftarrow CH-' name(4:5) ''',''Color'',''r'')']);
    try
    E_exp=pchip((EXP(:,2)),(abs(EXP(:,4))),displ);
    plot(displ, E_exp+0., '--k'); hold on;
    end
    % PLOT ALL FORCE
    figure(10);
    plot(displ, F_exp+0., '-k'); hold on;
    plot(displ, F_FEA+0., '-r'); hold on;
    eval(['text(EXP(end,2),EXP(end,3),''\leftarrow CH-' name(4:5) ''',''Color'',''k'')']); 
    eval(['text(displ(end),F_FEA(end),''\leftarrow CH-' name(4:5) ''',''Color'',''r'')']);
    
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
    %CHOOSE DURATION
    u_beg = EXP(end,2)*0.00;
    switch All_x_Fmax
    case 'All' 
        u_end = EXP(end,2);
    case 'Fmax'     
        [ic line]=max(EXP(:,3));
        u_end = max(EXP(line,2))  
    end
    figure(5);
    plot(EXP(:,2), EXP(:,3), '-g'); hold on;
    eval(['text(EXP(end,2),EXP(end,3),''\leftarrow SH-' name(4:5) ''',''Color'',''k'')']);  
    %---------------------------------------
    [B, t_sim, s] = readColData('temp.out',3,3);
    [kk, ia, ic]  = unique(s(:,1));
    ia=sort(ia);
    s = s(ia,:);
    displ = linspace(u_beg,u_end)';
    for j=1:75, coeff(j)=1;end
    for j=76:95, coeff(j)=1;end
    for j=96:100, coeff(j)=1;end
    %-----------------------------------
    % SHIFT SIM DISPL TO MATCH ELASTIC RANGE OF EXP
    eval(['displ=displ-deltaU.' name ';']);
    %-----------------------------------
    %Force pchip
	eval(['F_FEA=pchip((s(:,1)+deltaU.' name '),(-2./1000.*s(:,2)),displ)']);
    plot(displ, F_FEA+0., '-r'); hold on;
    eval(['text(displ(end),F_FEA(end),''\leftarrow SH-' name(4:5) ''',''Color'',''r'')']);
    F_exp=pchip((EXP(:,2)),(EXP(:,3)),displ);
    plot(displ, F_exp+0., '--k'); hold on;
    % PLOT ALL FORCE
    figure(10);
    plot(displ, F_exp+0., '-k'); hold on;
    plot(displ, F_FEA+0., '-r'); hold on;
    eval(['text(EXP(end,2),EXP(end,3),''\leftarrow SH-' name(4:5) ''',''Color'',''k'')']);
    eval(['text(displ(end),F_FEA(end),''\leftarrow SH-' name(4:5) ''',''Color'',''r'')']);
        
    err(i) = sum(coeff'.*(F_FEA-F_exp).^2);%+100.*sum(coeff'.*(E_FEA-E_exp).^2);
    eval(['disp('' MESSAGE: ' name ' POSTPROCESSING COMPLETED'')']);
    clear EXP F_FEA F_exp E_FEA E_exp displ s B t_sim ia ic j kk u_beg u_end coeff name; 
        
end
delete temp.out temp.odb
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Post Processing Completed')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%--------------------------------------------------------------
if ishandle(1)
figure(1);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
xlabel('Plastic Strain [-]','FontSize',12);
ylabel('True Stress [MPa]','FontSize',12);
yaxismin = 10000;
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
saveas(gcf,['_curve_' mat2str(Loop) '_UT.png'],'png');
% eval(['savefig(1,''_curve_' num2str(Loop) '_UT'')']);
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
saveas(gcf,['_curve_' mat2str(Loop) '_NT20.png'],'png');
% eval(['savefig(2,''_curve_' num2str(Loop) '_NT20'')']);
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
saveas(gcf,['_curve_' mat2str(Loop) '_NT6.png'],'png');
% eval(['savefig(3,''_curve_' num2str(Loop) '_NT6'')']);
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
% axis square;
box on
saveas(gcf,['_curve_' mat2str(Loop) '_CH.png'],'png');
% eval(['savefig(4,''_curve_' num2str(Loop) '_CH'')']);
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
saveas(gcf,['_curve_' mat2str(Loop) '_SH.png'],'png');
% eval(['savefig(5,''_curve_' num2str(Loop) '_SH'')']);
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
saveas(gcf,['_curve_' mat2str(Loop) '_AllFrac.png'],'png');
hold off;
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Image Save Complete')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%------------------------------------------------
% SUM ERROR AND FINISH POSTPROCESSING
%------------------------------------------------
result = sum(err);

namecsv = ['_material_iterations.csv'];
if Loop==0
fid=fopen(namecsv,'w');
fprintf(fid,'Loop,P12, P22, P44, G12, G22, G44, SValpha,\n');  
fclose(fid);
end

dlmwrite(namecsv,[Loop pv p12 p22 p44 g12 g22 g44 err result], 'precision', '%.6f', '-append')

if result < minres
    minres = result; 
    cd _minimum
    delete *.*
    cd ..  
    for i=1:nbExps  
    name=char(names(i))
    eval(['copyfile(''temp_' name '.odb'',''_minimum'')']);
    eval(['copyfile(''temp_' name '.out'',''_minimum'')']);
    eval(['copyfile(''temp_' name '.sta'',''_minimum'')']);
    
    if ishandle(1)
%     eval(['copyfile(''_curve_' num2str(Loop) '_UT.fig'',''_minimum'')']);
    eval(['copyfile(''_curve_' num2str(Loop) '_UT.png'',''_minimum'')']);
    end
    if ishandle(2)
%     eval(['copyfile(''_curve_' num2str(Loop) '_NT20.fig'',''_minimum'')']);
    eval(['copyfile(''_curve_' num2str(Loop) '_NT20.png'',''_minimum'')']);
    end
    if ishandle(3)
%     eval(['copyfile(''_curve_' num2str(Loop) '_NT6.fig'',''_minimum'')']);
    eval(['copyfile(''_curve_' num2str(Loop) '_NT6.png'',''_minimum'')']);
    end
    if ishandle(4)
%     eval(['copyfile(''_curve_' num2str(Loop) '_CH.fig'',''_minimum'')']);
    eval(['copyfile(''_curve_' num2str(Loop) '_CH.png'',''_minimum'')']);
    end
    if ishandle(5)
%     eval(['copyfile(''_curve_' num2str(Loop) '_SH.fig'',''_minimum'')']); 
    eval(['copyfile(''_curve_' num2str(Loop) '_SH.png'',''_minimum'')']); 
    end
    if ishandle(10)
    eval(['copyfile(''_curve_' num2str(Loop) '_AllFrac.png'',''_minimum'')']); 
    end
    copyfile('material.dat','_minimum')
    end
end  

%------------------------------------------------
% COPY OLD ITERATION RESULTS TO LAST FOLDER
%------------------------------------------------
    for i=1:nbExps  
    name=char(names(i))
    eval(['copyfile(''temp_' name '.odb'',''_previous'')']);
    eval(['copyfile(''temp_' name '.out'',''_previous'')']);
    eval(['copyfile(''temp_' name '.sta'',''_previous'')']);
    
    if ishandle(1)
%     eval(['copyfile(''_curve_' num2str(Loop) '_UT.fig'',''_minimum'')']);
    eval(['copyfile(''_curve_' num2str(Loop) '_UT.png'',''_previous'')']);
    end
    if ishandle(2)
%     eval(['copyfile(''_curve_' num2str(Loop) '_NT20.fig'',''_minimum'')']);
    eval(['copyfile(''_curve_' num2str(Loop) '_NT20.png'',''_previous'')']);
    end
    if ishandle(3)
%     eval(['copyfile(''_curve_' num2str(Loop) '_NT6.fig'',''_minimum'')']);
    eval(['copyfile(''_curve_' num2str(Loop) '_NT6.png'',''_previous'')']);
    end
    if ishandle(4)
%     eval(['copyfile(''_curve_' num2str(Loop) '_CH.fig'',''_minimum'')']);
    eval(['copyfile(''_curve_' num2str(Loop) '_CH.png'',''_previous'')']);
    end
    if ishandle(5)
%     eval(['copyfile(''_curve_' num2str(Loop) '_SH.fig'',''_minimum'')']); 
    eval(['copyfile(''_curve_' num2str(Loop) '_SH.png'',''_previous'')']); 
    end
    if ishandle(10)
    eval(['copyfile(''_curve_' num2str(Loop) '_AllFrac.png'',''_previous'')']); 
    end
    copyfile('material.dat','_previous')
    end

%------------------------------------------------
% MOVE ITERATION RESULTS TO ITERATIONS FOLDER
%------------------------------------------------
if ishandle(1)
% eval(['movefile(''_curve_' num2str(Loop) '_UT.fig'',''_iterations'')']);
% eval(['copyfile(''_curve_' num2str(Loop) '_UT.png'',''_iterations'')']);
eval(['movefile(''_curve_' num2str(Loop) '_UT.png'',''_overview/_iterations'')']);
end
if ishandle(2)
% eval(['movefile(''_curve_' num2str(Loop) '_NT20.fig'',''_iterations'')']);
eval(['movefile(''_curve_' num2str(Loop) '_NT20.png'',''_overview/_iterations'')']);
end
if ishandle(3)
% eval(['movefile(''_curve_' num2str(Loop) '_NT6.fig'',''_iterations'')']);
eval(['movefile(''_curve_' num2str(Loop) '_NT6.png'',''_overview/_iterations'')']);
end
if ishandle(4)
% eval(['movefile(''_curve_' num2str(Loop) '_CH.fig'',''_iterations'')']);
eval(['movefile(''_curve_' num2str(Loop) '_CH.png'',''_overview/_iterations'')']);
end
if ishandle(5)
% eval(['movefile(''_curve_' num2str(Loop) '_SH.fig'',''_iterations'')']); 
eval(['movefile(''_curve_' num2str(Loop) '_SH.png'',''_overview/_iterations'')']); 
end
if ishandle(10)
eval(['movefile(''_curve_' mat2str(Loop) '_AllFrac.png'',''_overview'')']);
end
copyfile('material.dat','material_old.dat')
hold off;
close all;

eval(['disp('' MESSAGE: LOOP #' num2str(Loop) ' COMPLETED'')'])
Loop = Loop+1;

clear a0 alpha ans E Etemp fid GG KK launchcount name NU spectemp Swift* Voce* 
clear SValpha p12 p22 p44 g12 g22 g44
exit
end
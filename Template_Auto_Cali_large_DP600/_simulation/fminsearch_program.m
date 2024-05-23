clear all; close all; format long; clc
global experiments deltaU nbExps YLD2000_pre Hill48_pre minres Loop;
delete screen_global.txt
diary screen_global.txt

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
    temp=csvread(files{i})
    tempU=temp(:,2)
    tempF=temp(:,3)
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
%-----------------------------------
YLD2000_pre = csvread('_YLD2000_pre.csv');
Hill48_pre = csvread('_Hill48naFr_pre.csv');
%-----------------------------------
% CREATE FOLDERS AND CLEANUP 
%-----------------------------------
mkdir _minimum
cd _minimum
delete *.*
cd ..

mkdir _overview
cd _overview
delete *.*
%
mkdir _iterations
cd _iterations
delete *.*
cd ..
%
cd ..
mkdir _previous
cd _previous
delete *.*
cd ..
%-----------------------------------
clear i ii listOf files name
%
Loop=0;
minres = 10E6;
%!del result.dat
options = optimset('Display','iter', 'TolFun',1e-4, 'TolX',1e-4);
% 
x0=[YLD2000_pre(7,7) 1.11 1.08 ];

% x0=[abs(Hill48_pre(5,4:6))];

% lb = [abs(Hill48_pre(5,4:6)*0.2) ];
% ub = [abs(Hill48_pre(5,4:6)*2.0) ];


lb = [ 0.0001 0.95 0.95];
ub = [ 0.9999 1.2 1.2];

% lb = [ 0.0001 ]; %0.5 0.3
% ub = [ 1.0000]; %3.0 0.7

% fminsearch
% [x,fval,exitflag,output] = fminsearch(@obj_function, x0, options) % Nelder_Mead direct search (simplex, no derivatives)

%fminsearchbnd
options = optimset('Display','iter','TolFun',1e-4, 'TolX',1e-4);
[x,fval,exitflag,output] = fminsearchbnd(@obj_function, x0, lb, ub, options)

%fmincon
% options = optimset('Display','iter', 'TolFun',1e-8, 'TolX',1e-8);
% [x,fval,exitflag,output] = fmincon(@obj_function, x0, [], [], [], [], lb, ub, [], options)
% 'DiffMinChange',5.e-2,'DiffMaxChange',1.e-1,








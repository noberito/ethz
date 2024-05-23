clear all; close all; format long; clc
%-----------------------------------
% LOAD AVAILABLE LOADPATHS 
%-----------------------------------
listOf = dir('*.csv');
ii=0;
for i=1:numel(listOf)
files{i}=listOf(i).name;
name=char(files(i));
name=name(1:end-4);
if any(regexp(name,'_path_'))==1
 eval(['loadpath.' name(7:end) '=csvread(files{i});']); 
 orientations{ii+1}=name(end-1:end);
 ii=ii+1;
 nbExps = ii;
end
end
orientations=sort(unique(orientations))';
names=fieldnames(loadpath);
clear i ii listOf files name

%-----------------------------------
% PREPARE ORIENTATIONS
%-----------------------------------
namecsv = ['_EMC-Parameters.csv'];
fid=fopen(namecsv,'w');
fprintf(fid,'Orientation,a,b,c,n,\n');  
fclose(fid);
    
for k=1:numel(orientations)
    clear ori count 
    clearvars -global filenames
    
    count=1;
    ori=char(orientations(k))
    
    global filenames

    for j=1:numel(fieldnames(loadpath))
        name=char(names(j))
        if strncmp(name(end-1:end),ori,2)==1
            eval(['filenames{count}=''' name ''';']);
            count=count+1;
        end
    end
    clear j

    % load  experimental data
    for i=1:length(filenames);
        clear BB CC;
        eval(['global ' char(filenames(i))]);
        eval(['BB =loadpath.' char(filenames(i)) '']);       
%         BB =csvread(strcat(char(filenames(i)),'.csv'));
        CC=EMC_prepare_data(BB);

        eval([char(filenames(i)) '=CC;']);
        
    end;   
    clear BB CC;
    
    % Prepare Minimization Parameters
    
    %[x,fval,exitflag,output] = fminsearch(@obj_function, x0, options) % Nelder_Mead direct search (simplex, no derivatives) 
    
    pv = [1.05 0.65 1.11];
    
    lb = [1.0 0.0 0.0];
    ub = [2.0 10. 20.];

    %fminsearchbnd
    options = optimset('Display','iter', 'TolFun',1e-7, 'TolX',1e-7)
    [xEMC,fval,exitflag,output] = fminsearchbnd(@EMC_obj_function, pv, lb, ub, options)
    
  
    %VISUALIZE 2D    
    %-------------------EMC------------------------------
    a=xEMC(1)
    b=xEMC(2)
    c=xEMC(3)/100.
    n=0.1 %xEMC(4)
    
    disp(['Orientation ' num2str(ori) '°'])
    disp(['a=' num2str(a) ' , b=' num2str(b) ' , c=' num2str(c) ' , n=' num2str(n) ''])
    
    dlmwrite(namecsv,[str2num(ori) a b c n], 'precision', '%.6f', '-append') 
    
    for j=1:length(filenames);
    clear eta D theta peeq X;
    X=eval(char(filenames(j)));
    imax=100;
    dpeeq = 5./imax;
    peeq(1) = 0.;
    D(1)=0.;
    for i=2:imax
        peeq(i)  = peeq(i-1) + dpeeq;
        eta(i)   = interp1(X(:,3), X(:,1),peeq(i),'linear');
        theta(i) = interp1(X(:,3), X(:,2),peeq(i),'linear');
        epsf     = EMC(eta(i), theta(i),a,b,c,n);
        D(i)     = D(i-1) + dpeeq/epsf;
    end;
    %
        eps_frac  = interp1(D, peeq, 1.0,'linear')
        eta_frac  = interp1(D, eta, 1.0,'linear');
        theta_frac  = interp1(D, theta, 1.0,'linear');
    ll = length(X)-1;
    % PLOT INDIVIDUAL ORIENTATIONS
    figure(str2num(ori)+101); %TRIAX
    plot(X(1:ll,1), X(1:ll,3),'k-'); hold on;
    plot(eta_frac, eps_frac,'ok', 'MarkerSize', 5, 'MarkerEdgeColor','k','MarkerFaceColor','b'); hold on;
    figure(str2num(ori)+102); %LODE
    plot(X(1:ll,2), X(1:ll,3),'k-'); hold on;
    plot(theta_frac, eps_frac,'ok', 'MarkerSize', 5, 'MarkerEdgeColor','k','MarkerFaceColor','b'); hold on;    
    % PLOT ALL ORIENTATIONS
    figure(31); %TRIAX
    plot(X(1:ll,1), X(1:ll,3),'k-'); hold on;
    plot(eta_frac, eps_frac,'ok', 'MarkerSize', 5, 'MarkerEdgeColor','k','MarkerFaceColor','b'); hold on;
    figure(32); %LODE
    plot(X(1:ll,2), X(1:ll,3),'k-'); hold on;
    plot(theta_frac, eps_frac,'ok', 'MarkerSize', 5, 'MarkerEdgeColor','k','MarkerFaceColor','b'); hold on;
    xmax(k,j)=max(X(1:ll,1));
    ymax(k,j)=max(X(1:ll,3));
    end;
    %--------------------------------------------------------------------------
    % plane stress visualization
    clear eta xi theta epsf
    eta= [-1/3:1/300:2/3]; % eta
    xi=-27/2*eta.*(eta.^2-1/3);
    theta=1-2/pi*acos(xi);
    for i=1:length(eta);
        epsf(i)=EMC(eta(i), theta(i), a, b, c,n);
    end;
    % PLOT INDIVIDUAL ORIENTATIONS
    figure(str2num(ori)+101); %TRIAX
    plot(eta, epsf,'b-'); hold on;
    eval(['text(eta(end),epsf(end),''\leftarrow ' ori ''',''Color'',''k'')']);     
    figure(str2num(ori)+102); %LODE 
    plot(theta, epsf,'b-'); hold on;
    eval(['text(theta(end-1),epsf(end-1),''\leftarrow ' ori ''',''Color'',''k'')']);    
    % PLOT ALL ORIENTATIONS
    figure(31); %TRIAX
    plot(eta, epsf,'b-'); hold on;
    eval(['text(eta(end),epsf(end),''\leftarrow ' ori ''',''Color'',''k'')']); 
    figure(32); %LODE
    plot(theta, epsf,'b-'); hold on;
    eval(['text(theta(end-1),epsf(end-1),''\leftarrow ' ori ''',''Color'',''k'')']);   
end   
xmax=max(max(xmax));
ymax=max(max(ymax));
%--------------------------------------------------------------
for lll=[31 101 191]
   try
    figure(lll);
    set(gca,'XMinorTick','on');
    set(gca,'YMinorTick','on');
    set(gca,'TickLength',[0.01 0.01]);
    set(gca,'FontSize',12);
    set(gca,'LineWidth',1);
    set(gca,'XTick',[-0.33,-.2,-.1,0,0.1,0.2,0.33,0.45,0.57,0.67,0.8,0.9,1]);
    xlabel('Triaxiality [-]','FontSize',12);
    ylabel('Eq. Plastic Strain [-]','FontSize',12);
    axis([-0.1 1.1*xmax 0.0 1.1*ymax])
    axis square;
    box on
        if lll == 31 
        saveas(gcf,['_final_Triax_PEEQ_HC.png'],'png');
        savefig(lll,'_final_Triax_PEEQ_HC');
        elseif lll == 101 
        saveas(gcf,['_final_Triax_PEEQ_HC_00.png'],'png');
        savefig(lll,'_final_Triax_PEEQ_HC_00');
        elseif lll == 191 
        saveas(gcf,['_final_Triax_PEEQ_HC_90.png'],'png');
        savefig(lll,'_final_Triax_PEEQ_HC_90');
        end
    hold off;
    end
end
%--------------------------------------------------------------
for lll=[32 102 192]
    try figure(lll);
    set(gca,'XMinorTick','on');
    set(gca,'YMinorTick','on');
    set(gca,'TickLength',[0.01 0.01]);
    set(gca,'FontSize',12);
    set(gca,'LineWidth',1);
    xlabel('Lode Angle Parameter [-]','FontSize',12);
    ylabel('Eq. Plastic Strain [-]','FontSize',12);
    axis([-1.0 1.0 0.0 1.1*ymax])
    axis square;
    box on
        if lll == 32 
        saveas(gcf,['_final_Lode_PEEQ_HC.png'],'png');
        savefig(lll,'_final_Lode_PEEQ_HC');
        elseif lll == 102 
        saveas(gcf,['_final_Lode_PEEQ_HC_00.png'],'png');
        savefig(lll,'_final_Lode_PEEQ_HC_00');
        elseif lll == 192 
        saveas(gcf,['_final_Lode_PEEQ_HC_90.png'],'png');
        savefig(lll,'_final_Lode_PEEQ_HC_90');
        end
    hold off;
    hold off;
    end
end








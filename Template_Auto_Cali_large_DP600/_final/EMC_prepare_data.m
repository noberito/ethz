function result = prepare_data(BB);
% input: BB, col #1 = eta
%            col #2 = theta_bar
%            col #3 = peeq
%
% output: same matrix with 
%           - all rows of zero straion increments deleted
%           - an additional line for a strain of 5
%
    BB(1,4)=1;
    for j=2:length(BB)
        BB(j,4)=BB(j,3)-BB(j-1,3);
        if (BB(j,4)<0) 
            BB(j,4)=0.; 
        end;
    end;
    BB
    % 2. delete rows of zero strain increment
    bb = any(BB(:,4)~=0,2);
    BB = BB(bb,:);
    [C,IA,IC] = unique(BB(:,3),'rows','stable');
    BB = BB(IA,:);
    clear C IA IC
    
    peeq_temp = linspace(BB(1,3),BB(end,3),1000);
    BB1_temp = interp1(BB(:,3),BB(:,1),peeq_temp);
    BB2_temp = interp1(BB(:,3),BB(:,2),peeq_temp);
    clear BB
    BB(:,1)  = BB1_temp;
    BB(:,2)  = BB2_temp;
    BB(:,3)  = peeq_temp;
    clear BB1_temp BB2_temp peeq_temp 
    
    % 3. add line
    bl = length(BB);
    % extrapolation for stress triaxiality
    % x = peeq; y = eta;
        x_b = BB(bl,3);
        x_a = 0.8*x_b;
        y_b = BB(bl,1);
        y_a = interp1(BB(:,3), BB(:,1),x_a,'linear');
    %
        slope = (y_b-y_a)/(x_b-x_a);
        x_c = 5.;
        y_c = y_b + slope*(x_c-x_b);
    %
        BB(bl+1, 1)=y_c; 
        BB(bl+1, 3)=x_c; 
%-------------------------------------------------
    % extrapolation for Lode parameter
    % x = peeq; y = theta;
        x_b = BB(bl,3);
        x_a = 0.8*x_b;
        y_b = BB(bl,2);
        y_a = interp1(BB(:,3), BB(:,2),x_a,'linear');
    %
        slope = (y_b-y_a)/(x_b-x_a);
        x_c = 5.;
        y_c = y_b + slope*(x_c-x_b);
    %
        if(y_c>1) y_c=1; end;
        if(y_c<-1) y_c=-1; end;
        BB(bl+1, 2)=y_c; 
%-------------------------------------------------
result=BB(:,1:3);
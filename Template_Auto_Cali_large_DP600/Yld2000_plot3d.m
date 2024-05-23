function Sigma= Yld2000_plot3d(X,Y,Z,alpha,a)

[L1,L2] = Yld2000_getMatrices(alpha); 

Sigma = nan(size(X)); 


for i=1:size(X,1)
    
    for j=1:size(X,2)
  
        Stress_State = [X(i,j),Y(i,j),Z(i,j)]';
        sigma0 = 1;
                
        options = optimset('Display', 'off');

        Sigma(i,j) = fsolve(@(sigma0) plot_err(sigma0,Stress_State, a,L1,L2), sigma0, options);
%         Sigma(i,j) = fsolve(@(sigma0) plot_err(sigma0,Stress_State, alpha, a), sigma0, options);


    end
    
end

end

% function err = plot_err(sigma,Stress_State, alpha, a)
%     sig_vec = sigma*Stress_State;
%     err = yld2000_fun_old(sig_vec, alpha, a);
% end


function err = plot_err(sigma,Stress_State, a,L1,L2)
    sig_vec = sigma*Stress_State;
    err = yld2000_fun(sig_vec, a,L1,L2);
end

function [L1,L2] = Yld2000_getMatrices(alpha)
    M1 = [2/3 0 0; -1/3 0 0; 0 -1/3 0; 0 2/3 0; 0 0 1];
    M2 = 1/9*[-2 2 8 -2 0; 1 -4 -4 4 0; 4 -4 -4 1 0; -2 8 2 -2 0; 0 0 0 0 9];

    l1 = M1 * [alpha(1); alpha(2); alpha(7)];
    l2 = M2 * [alpha(3); alpha(4); alpha(5); alpha(6); alpha(8)];

    L1 = [l1(1) l1(2) 0; l1(3) l1(4) 0; 0 0 l1(5)];
    L2 = [l2(1) l2(2) 0; l2(3) l2(4) 0; 0 0 l2(5)];
end



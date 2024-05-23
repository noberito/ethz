function result = obj_function_MC(pv);
%
global filenames sig eps refx curpara;
for i=1:length(filenames)
    eval(['global ' char(filenames(i))]);
end;
% EMC:
 a=pv(1); 
 b=pv(2);
 c=pv(3)/100; 
 n=0.1;
for j=1:length(filenames);
    clear eta D theta peeq X;
    X=eval(char(filenames(j)));
    % EMC(0.33, 1, 1.47 1020.8, 0.008)
    imax=100;
    dpeeq = 2./imax;
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
    if(max(D)<1)
        eps_frac=5.;
    else
        eps_frac  = interp1(D, peeq, 1.0,'linear');
    end;
    eps_exp = X(length(X)-1,3);
    err(j) = abs(eps_frac/eps_exp-1.);
end;
result = max(abs(err));
%result = sum(err)/length(err);
out2 = [a b c result];
fprintf('%8.4f, %8.4f,%8.4f,%8.4f, \n', out2)
out = [a b c err];
curpara=[a b c];
%
fid=fopen('result.dat','a');
fprintf(fid, '\n %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f', out);
fclose(fid);
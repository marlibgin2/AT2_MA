function [X,fval,exitflag,output] = Optim_maxmed_OCT(X0,Rin)
pp = gcp;

global Rin 
if nargin < 1
    X0(1) = atgetfieldvalues(Rin(oxxoi(1)),'PolynomB',{1,4}); 
    X0(2) = atgetfieldvalues(Rin(oxyoi(1)),'PolynomB',{1,4}); 
    X0(3) = atgetfieldvalues(Rin(oyyoi(1)),'PolynomB',{1,4}); 
end

delta = [500, 1000 ,500];
options = optimset('Display','iter','MaxIter',150,'MaxFunEvals',250,'TolFun',1e-5,'TolX',1e-5,'PlotFcns',@optimplotfval);

[X,fval,exitflag,output] = fmincon(@fun_maxmed_OCT,X0, [], [], [], [], X0-delta, X0+delta, [], options);

delete(pp)
end
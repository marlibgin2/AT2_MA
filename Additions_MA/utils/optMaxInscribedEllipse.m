function [X,fval,exitflag,output] = optMaxInscribedEllipse(X0)
options = optimset('Display','iter','MaxIter',300,'MaxFunEvals',300,'TolFun',1e-4,'TolX',1e-4,'PlotFcns',@optimplotfval);


delta = [10,10,10];
[X,fval,exitflag,output] = fminsearchbnd(@funMaxInscribedEllipse,X0, X0-delta, X0+delta, options);
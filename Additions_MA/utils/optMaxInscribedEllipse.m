function [X,fval,exitflag,output] = optMaxInscribedEllipse(X0)
options = optimset('Display','iter','MaxIter',400,'MaxFunEvals',400,'TolFun',1e-6,'TolX',1e-6)%,'PlotFcns',@optimplotfval);


delta = [10,10,10];
%[X,fval,exitflag,output] = fminsearchbnd(@funMaxInscribedEllipse,X0, X0-delta, X0+delta, options);
[X,fval,exitflag,output] = fminsearchcon(@funMaxInscribedEllipse,X0, X0-delta, X0+delta, [], [], @nonlcon, options);
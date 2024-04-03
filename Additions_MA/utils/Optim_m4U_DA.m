% ---------------------------------------------------------
% optimisation of DA of a given RING 
% ---------------------------------------------------------
function [daPAR,fval,exitflag,output] = Optim_m4U_DA(daPAR0)
%options = optimset('Display','iter','MaxIter',500,'MaxFunEvals',500,'TolFun',1e-3,'TolX',1e-3,'PlotFcns',@optimplotfval);
options = optimset('Display','iter','MaxIter',50,'MaxFunEvals',100,'TolFun',1e-3,'TolX',1e-3,'PlotFcns',@optimplotfval);

% nux, nuy, oxxo, oxyo oyyo
dPAR = [0.23, 0.23, 2000, 2000, 2000];
%delta = [0.1, 0.1, 0.05, 0.05, 0.005,   0.005, 0.001, 0.002, 0.002, 0.002, 0.005 ]/10;
ulPAR = [56.35 16.35  2000  2000  2000];
llPAR = [56.05 16.05 -2000 -2000 -2000];
[daPAR,fval,exitflag,output] = fminsearchbnd(@fun_m4U_DA_AT2,daPAR0, llPAR, ulPAR, options);


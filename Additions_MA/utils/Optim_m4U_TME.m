% ---------------------------------------------------------
% optimisation of a mTME cell
% ---------------------------------------------------------
function [X,fval,exitflag,output] = Optim_m4U_TME(X0)
%options = optimset('Display','iter','MaxIter',500,'MaxFunEvals',500,'TolFun',1e-3,'TolX',1e-3,'PlotFcns',@optimplotfval);
options = optimset('Display','iter','MaxIter',10,'MaxFunEvals',20,'TolFun',1e-3,'TolX',1e-3,'PlotFcns',@optimplotfval);


delta = [0.1, 0.1];
delta(1) = .05*8; delta(2) = .05*8; delta(3) = .005; 
%delta(3) = 0.2; delta(4) = 0.2; delta(5) = 0.2; delta(6) = 0.4; 
%delta(1:6) = delta(1:6)/4;

%delta = [0.1, 0.1, 0.05, 0.05, 0.005,   0.005, 0.001, 0.002, 0.002, 0.002, 0.005 ]/10;

[X,fval,exitflag,output] = fminsearchbnd(@fun_m4U_TME_match_AT2,X0, X0-delta, X0+delta, options);


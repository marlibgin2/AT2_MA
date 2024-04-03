% ---------------------------------------------------------
% optimisation of a mTME cell
% ---------------------------------------------------------
function [X,fval,exitflag,output] = Optim_numax(X0)
global rrr
%rrr   = numax_7BA_25_2_4; 
load('numax_7BA_25_2_4_M1.mat','RING','RING1');rrr=RING1; clear RING RING1;

if nargin < 1
%X0 = [33.511620518476079 -18.537623834316843 -5.947869392024406 -315.3 241.400 -94.22];
%X0 = [0.659875782646984  -0.647079937939663   0.094501286907741 -2.932012994574695  -2.270273823218202  -0.921669116538125]*100

hsi = findcells(rrr,'FamName','HS');
hsfi = findcells(rrr,'FamName','HSF');
hsdi = findcells(rrr,'FamName','HSD');

oxxoi  = findcells(rrr,'FamName','oxxo');
oxyoi  = findcells(rrr,'FamName','oxyo');
oyyoi  = findcells(rrr,'FamName','oyyo');

X0(1) = atgetfieldvalues(rrr(hsi(1)),'PolynomB',{1,3}); 
X0(2) = atgetfieldvalues(rrr(hsfi(1)),'PolynomB',{1,3}); 
X0(3) = atgetfieldvalues(rrr(hsdi(1)),'PolynomB',{1,3});  

X0(4) =atgetfieldvalues(rrr(oxxoi(1)),'PolynomB',{1,4}); 
X0(5) = atgetfieldvalues(rrr(oxyoi(1)),'PolynomB',{1,4}); 
X0(6) = atgetfieldvalues(rrr(oyyoi(1)),'PolynomB',{1,4}); 
end

options = optimset('Display','iter','MaxIter',60,'MaxFunEvals',100,'TolFun',1e-3,'TolX',1e-3,'PlotFcns',@optimplotfval);

%        HS, HSF, HSD, oxxxo, oxyo, oyyo 
delta = [200,  200,  200,   400,  400,  400];


[X,fval,exitflag,output] = fminsearchbnd(@fun_ax_AT2,X0, X0-delta, X0+delta, options);

hsi = findcells(rrr,'FamName','HS');
hsfi = findcells(rrr,'FamName','HSF');
hsdi = findcells(rrr,'FamName','HSD');

oxxoi  = findcells(rrr,'FamName','oxxo');
oxyoi  = findcells(rrr,'FamName','oxyo');
oyyoi  = findcells(rrr,'FamName','oyyo');

for i = 1:length(hsi)
    rrr{hsi(i)}.PolynomB(3) = X(1);
end
for i = 1:length(hsfi)
    rrr{hsfi(i)}.PolynomB(3) = X(2);
end
for i = 1:length(hsdi)
    rrr{hsdi(i)}.PolynomB(3) = X(3);
end

for i = 1:length(oxxoi)
    rrr{oxxoi(i)}.PolynomB(4) = X(4);
end
for i = 1:length(oxyoi)
    rrr{oxyoi(i)}.PolynomB(4) = X(5);
end
for i = 1:length(oyyoi)
    rrr{oyyoi(i)}.PolynomB(4) = X(6);
end

tic; [xf_da, yf_da, DAf] = moga_DA(rrr,15,40,500,7e-3,3e-3,0.0);toc % the triangular DA! 
disp(['DA = ' num2str(DAf)])
figure(112); hold on; plot(xf_da, yf_da, '-o','linewidth',3); 



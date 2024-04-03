% ---------------------------------------------------------
% optimisation of a mTME cell
% ---------------------------------------------------------
function [X,fval,exitflag,output] = Optim_medmax(X0)
global rrr
%rrr   = numax_7BA_25_2_4; 
%load('numax_7BA_25_2_4_M1.mat','RING','RING1');rrr=RING1; clear RING RING1;

%load('medmax_7BA_1_1_1_AT2_simple.mat','RING','RING1');rrr=RING1; clear RING RING1;
load('medmax_7BA_2_1_1_AT2_simple.mat','RING','RING1');rrr=RING1; clear RING RING1;
%load('medmax_7BA_2_5_1_AT2_simple.mat','RING','RING1');rrr=RING1; clear RING RING1;


%qfi    = findcells(rrr,'FamName','QF');
qfii    = findcells(rrr,'FamName','QFI');
qfoi    = findcells(rrr,'FamName','QFO');

qfmi   = findcells(rrr,'FamName','QFM');
qfendi = findcells(rrr,'FamName','QFEND');
qdendi = findcells(rrr,'FamName','QDEND');

sdi     = findcells(rrr,'FamName','SD');
sdendi  = findcells(rrr,'FamName','SDEND');
sfmi    = findcells(rrr,'FamName','SFM');
sfoi    = findcells(rrr,'FamName','SFO');
sfii    = findcells(rrr,'FamName','SFI');

dipi   = findcells(rrr,'FamName','DIP');
dipmi  = findcells(rrr,'FamName','DIPm');


for i=1:length(dipi)
    wdip(i) = rrr{dipi(i)}.PolynomB(2)/rrr{dipi(7)}.PolynomB(2);
end
for i=1:length(dipmi)
    wdipm(i) = rrr{dipmi(i)}.PolynomB(2)/rrr{dipmi(7)}.PolynomB(2);
end

if nargin < 1
X0(1) = atgetfieldvalues(rrr(qfii(1)),'PolynomB',{1,2}); 
X0(2) = atgetfieldvalues(rrr(qfoi(1)),'PolynomB',{1,2}); 
X0(3) = atgetfieldvalues(rrr(qfmi(1)),'PolynomB',{1,2}); 
X0(4) = atgetfieldvalues(rrr(qfendi(1)),'PolynomB',{1,2});  
X0(5) = atgetfieldvalues(rrr(qdendi(1)),'PolynomB',{1,2});  
X0(6) = atgetfieldvalues(rrr(dipi(7)),'PolynomB',{1,2}); 
X0(7) = atgetfieldvalues(rrr(dipmi(7)),'PolynomB',{1,2}); 
X0(8) = atgetfieldvalues(rrr(sdi(1)),'PolynomB',{1,3}); 
X0(9) = atgetfieldvalues(rrr(sdendi(1)),'PolynomB',{1,3}); 
X0(10) = atgetfieldvalues(rrr(sfmi(1)),'PolynomB',{1,3}); 
X0(11) = atgetfieldvalues(rrr(sfoi(1)),'PolynomB',{1,3}); 
X0(12) = atgetfieldvalues(rrr(sfii(1)),'PolynomB',{1,3}); 
end

options = optimset('Display','iter','MaxIter',150,'MaxFunEvals',250,'TolFun',1e-5,'TolX',1e-5,'PlotFcns',@optimplotfval);
% options = optimset('Display','iter','MaxIter',70,'MaxFunEvals',140,'TolFun',1e-4,'TolX',1e-4,'UseParallel',true);

%%% ref: delta = [0.01,  0.01, 0.01,  0.01,  0.01,  0.002,  0.002]*10;
%%% delta = [0.05,  0.05, 0.05,  0.04,  0.03,  0.004,  0.005]*10;
%%% delta = [0.05,  0.05, 0.05,  0.04,  0.03,  0.004,  0.01 2 2 2 2 2]*10;
delta = [0.05,  0.05, 0.05,  0.05,  0.05,  0.01,  0.01, 0, 0, 0, 0, 0]*50;


%%% [X,fval,exitflag,output] = fminsearchbnd(@fun_medmax_match_AT2,X0, X0-delta, X0+delta, options);
[X,fval,exitflag,output] = fmincon(@fun_medmax_match_AT2,X0, [], [], [], [], X0-delta, X0+delta, [], options);

if 1==0
    
for i = 1:length(qfi)
    rrr{qfi(i)}.PolynomB(2) = X(1);
end
for i = 1:length(qdi)
    rrr{qdi(i)}.PolynomB(2) = X(2);
end
for i = 1:length(qfendi)
    rrr{qfendi(i)}.PolynomB(2) = X(3);
end

for i = 1:length(qdendi)
    rrr{qdendi(i)}.PolynomB(2) = X(4);
end
for i = 1:length(dipi)
    rrr{dipi(i)}.PolynomB(2) = X(5).*wdipi(i);
end
for i = 1:length(dipmi)
    rrr{dipmi(i)}.PolynomB(2) = X(6).*wdipmi(i);
end


end
% ---------------------------------------------------------
% optimisation of a mTME cell
% ---------------------------------------------------------
function [X,fval,exitflag,output] = MOGA_medmax(X0)%, besti)
%global rrr

%rrr   = numax_7BA_25_2_4; 
%load('numax_7BA_25_2_4_M1.mat','RING','RING1');rrr=RING1; clear RING RING1;

%load('medmax_7BA_1_1_1_AT2_simple.mat','RING','RING1');rrr=RING1; clear RING RING1;
load('medmax_7BA_2_1_1_AT2_simple.mat','RING','RING1');rrr=RING1; clear RING RING1;
%load('medmax_7BA_2_5_1_AT2_simple.mat','RING','RING1');rrr=RING1; clear RING RING1;
Crrr = parallel.pool.Constant(rrr);


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

oxxoi = findcells(rrr,'FamName','oxxo');
oxyoi = findcells(rrr,'FamName','oxyo');
oyyoi = findcells(rrr,'FamName','oyyo');
VARi = {oxxoi; oxyoi; oyyoi};

% zero all the octupoles )in force from 24-11-2023)
noOCT=1;
if noOCT==1
for j=1:3
    for i=1:length(VARi{j})
        rrr{VARi{j}(i)}.PolynomB(4) = 0;
    end
end
end

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

% options = optimset('Display','iter','MaxIter',150,'MaxFunEvals',250,'TolFun',1e-5,'TolX',1e-5,'PlotFcns',@optimplotfval);
% options = optimset('Display','iter','MaxIter',70,'MaxFunEvals',140,'TolFun',1e-4,'TolX',1e-4,'UseParallel',true);

% % % options = optimoptions('gamultiobj','PlotFcn',@gaplotpareto,'PopulationSize',5000,'UseParallel',true,'UseVectorized', false,...
% % %                        'FunctionTolerance',1e-3,'MaxStallGenerations',5,'MaxGenerations',5, ...
% % %                        'ConstraintTolerance',1e-1);
%delta = [0.05,  0.05, 0.05,  0.05,  0.05,  0.01,  0.01, 0, 0, 0, 0, 0]*10;
%%%%%%% delta  = [0.05,  0.05, 0.05,  0.05,  0.05,  0.01,  0.01, 0.2, 0.2, 0.2, 0.2, 0.2];
%%%%%%% delta  = [0.44,  0.44, 0.38,  0.38,  0.25,  0.09,  0.08, 11, 16, 17, 17, 20]*1; % 10% variation
%%%%%%%delta  = [0.44,  0.44, 0.38,  0.38,  0.25,  0.09,  0.08, 11, 16, 17, 17, 20]*2; % 20% variation

%%%%% delta  = [0.44,  0.44, 0.38,  0.38,  0.25,  0.09,  0.08, 11*2, 16*2, 17*2, 17*2, 20*2]*1; % 10% variation for quad/Dipgrad, 20% variation for sext
%%%%% delta  = [0.44,  0.44, 0.38,  0.38,  0.25,  0.09,  0.08, 11*3, 16*3, 17*3, 17*3, 20*3]*1; % 10% variation for quad/Dipgrad, 30% variation for sext

%% trying to buzz around 155 pm.rad / I reduce the quad span and increase the sextupole changes
delta  = [0.44/2,  0.44/2, 0.38/2,  0.38/2,  0.25/2,  0.09/2,  0.08/2, 11*6, 16*6, 17*6, 17*6, 20*6]*1; % 5% variation for quad/Dipgrad, 60% variation for sext

sizPop = 10000; 
iniPop = 2*(rand(sizPop,12)-0.5).*delta/1+X0;
% iniPop = [iniPop; X0]; 

%sizIni = 100;
%iniPop = 2*(rand(sizIni,12)-0.5).*delta/100+X0;
%iniPop =  repmat(X0,sizPop,1);


% iniPop = 2*(rand(5000-size(X0,1),12)-0.5).*delta+X0(besti,:);
% iniPop = [iniPop; X0]; 

%iniPop = X0; 
controlgcp = 1; 
usepar = true;
if controlgcp==1
   c = parcluster;
   % next line is optional at MAX IV
   c.AdditionalProperties.AccountName = 'any-virtual-account-name';
   % 6 hour walltime
   c.AdditionalProperties.WallTime = '70:00:00';
   % hyperthreading enabled
   c.saveProfile
   c.NumThreads = 2;
   parpool('aurora R2022a',48) %48 %96
% % % %    c.NumThreads = 1;
% % % %    parpool('local',12) %48
   pp = gcp; 
   usepar = true;
end

options = optimoptions('gamultiobj','PlotFcn',@gaplotpareto,'PopulationSize',sizPop,'InitialPopulationMatrix',iniPop,...
                       'UseParallel',usepar,'UseVectorized', false,...
                       'display', 'iter', 'OutputFcns', @MOGAoutputfcn, ...
                       'FunctionTolerance',1e-3,'MaxStallGenerations',5,'MaxGenerations',20, ...
                       'ConstraintTolerance',1e-1);
%%% [X,fval,exitflag,output] = fminsearchbnd(@fun_medmax_match_AT2,X0, X0-delta, X0+delta, options);X
%%% [X,fval,exitflag,output] = fmincon(@fun_medmax_match_AT2,X0, [], [], [], [], X0-delta, X0+delta, [], options);

%fitnessfcn = @(X)[obj_emix(X, rrr), obj_rdt(X, rrr) ];
%%%%% fitnessfcn = @(X)[obj_emix(X, rrr), obj_modelDA(X, rrr) ];
fitnessfcn = @(X)OBJfcn(X, rrr); % object function
mycon      = @(X)CONfcn(X, rrr); % constraint function
nvars      = length(X0);
[X,fval,exitflag,output] = gamultiobj(fitnessfcn,nvars,[],[],[],[],X0-delta, X0+delta, mycon, options);

adesso = datestr(now,'ddmmyyyy_HHMMSS');
save(['MOGA_' adesso '.mat'],'X','fval','exitflag','output')

if controlgcp==1
    delete(pp)
end

end

% % % function obj_emix(X)
% % % global rrr
% % % 
% % % qfii    = findcells(rrr,'FamName','QFI');
% % % qfoi    = findcells(rrr,'FamName','QFO');
% % % 
% % % qfmi   = findcells(rrr,'FamName','QFM');
% % % qfendi = findcells(rrr,'FamName','QFEND');
% % % qdendi = findcells(rrr,'FamName','QDEND');
% % % 
% % % sdi     = findcells(rrr,'FamName','SD');
% % % sdendi  = findcells(rrr,'FamName','SDEND');
% % % sfmi    = findcells(rrr,'FamName','SFM');
% % % sfoi    = findcells(rrr,'FamName','SFO');
% % % sfii    = findcells(rrr,'FamName','SFI');
% % % 
% % % dipi   = findcells(rrr,'FamName','DIP');
% % % 
% % % dipmi  = findcells(rrr,'FamName','DIPm');
% % % 
% % % 
% % % for i=1:length(dipi)
% % %     wdip(i) = rrr{dipi(i)}.PolynomB(2)/rrr{dipi(7)}.PolynomB(2);
% % % end
% % % for i=1:length(dipmi)
% % %     wdipm(i) = rrr{dipmi(i)}.PolynomB(2)/rrr{dipmi(7)}.PolynomB(2);
% % % end
% % % 
% % % 
% % % 
% % % 
% % % end
% % % 
% % % function obj_rdt(X)
% % % global rrr
% % % 
% % % end
% % % 

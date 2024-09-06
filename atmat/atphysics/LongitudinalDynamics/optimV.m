function Voptimum = optimV(varargin)
% Calculates the optimum rf volatge needed to get maximum Touscheck
% lifetime 
% 
% Note (2024/09/07) : a small discrepancy (~10%) has been found between the model
% used in this function and the results of straightforward 6D tracking,
% under investigatio
%
%% Mandatory Inputs
% RING: Lattice array structure for whole ring
%
%% Outputs
% Voptimum : RF voltage [V]
%
%% Usage examples
% Voptimum=optimV(RING);

%% History
% S Jena 2024/08/19 : solve for Vopt for maximum Touscheck lifetime
% S Jena 2024/08/22 : changes made to use calc_TLT rather than calcTLT_raw

%% Input argument parsing
 [RING, TLoptions,MAoptions]= getargs(varargin,[],[],[]);

if (isempty(MAoptions))
    MAoptions.nperiods=20;
    MAoptions.lmafams='all';
    MAoptions.stepfam=1;
    MAoptions.deltalimit=0.3;
    MAoptions.initcoord=[0 0 0 0 0 nan]';
    MAoptions.delta=0.01;
    MAoptions.deltastepsize=0.1;
    MAoptions.splits=10;
    MAoptions.split_step_divisor=2;
    MAoptions.nturns=nan;
    MAoptions.S0max=528/20;
    MAoptions.S0min=0.0;
end

nperiods           = getoption(varargin,'nperiods',20);
lmafams            = getoption(varargin,'lmafams',MAoptions.lmafams);
stepfam            = getoption(varargin,'stepfam',MAoptions.stepfam);
deltalimit         = getoption(varargin,'deltalimit',MAoptions.deltalimit);
initcoord          = getoption(varargin,'initcoord',MAoptions.initcoord);
delta              = getoption(varargin,'delta',MAoptions.delta);
deltastepsize      = getoption(varargin,'deltastepsize',MAoptions.deltastepsize);
splits             = getoption(varargin,'splits',MAoptions.splits);
split_step_divisor = getoption(varargin,'split_step_divisor',MAoptions.split_step_divisor);
nturns             = getoption(varargin,'nturns',MAoptions.nturns);
S0max              = getoption(varargin,'S0max', MAoptions.S0max);
S0min              = getoption(varargin,'S0min', MAoptions.S0min);

MAoptions.nperiods   = nperiods;
MAoptions.lmafams    = lmafams;
MAoptions.stepfam    = stepfam;
MAoptions.deltalimit = deltalimit;
MAoptions.initcoord  = initcoord;
MAoptions.delta      = delta;
MAoptions.deltastepsize =  deltastepsize;
MAoptions.splits     = splits;
MAoptions.split_step_divisor = split_step_divisor;
MAoptions.nturns             = nturns;
MAoptions.S0max              = S0max;
MAoptions.S0min              = S0min;
% 
if isempty(TLoptions)
    TLoptions.Ib = 0.5/176;
    TLoptions.integrationmethod = 'integral';
    TLoptions.abstol = 1.0e-16;
    TLoptions.reltol = 1.0e-16;
    TLoptions.nperiods = 20;
    TLoptions.LMAperiods = 1;
    TLoptions.kcoupl = 'auto';
    TLoptions.emity  = 8.0E-12;
    TLoptions.OptimOvervoltage=100;
end
% verboselevel  = getoption(varargin,'verbose',0);
Ib = getoption(varargin,'Ib',TLoptions.Ib);
integrationmethod = getoption(varargin,'integrationmethod',TLoptions.integrationmethod);
abstol = getoption(varargin, 'AbsTol', TLoptions.abstol);
reltol = getoption(varargin, 'Relol', TLoptions.reltol);
nperiods = getoption(varargin,'nperiods',TLoptions.nperiods);
LMAperiods = getoption(varargin,'LMAperiods',TLoptions.LMAperiods);
kcoupl = getoption(varargin,'kcoupl',TLoptions.kcoupl);
emity  = getoption(varargin,'emity',TLoptions.emity);
OptimOvervoltage = getoption(varargin,'OptimOvervoltage',TLoptions.OptimOvervoltage);
harm = getoption(varargin,'harm',176);
% LMA    = getoption(varargin,'LMA', struct());
   
TLoptions.Ib=Ib;
TLoptions.integrationmethod = integrationmethod;
TLoptions.abstol = abstol;
TLoptions.reltol = reltol;
TLoptions.nperiods = nperiods;
TLoptions.LMAperiods = LMAperiods;
TLoptions.kcoupl = kcoupl;
TLoptions.emity  = emity;
TLoptions.OptimOvervoltage=OptimOvervoltage;

%% LMA calculation at high RF voltage 
U0=atgetU0(RING);
V0=TLoptions.OptimOvervoltage*U0;
alphac=mcf(RING);
energy=atenergy(RING);
cavpts  = find(atgetcells(RING, 'FamName', 'CAV'));
if (isempty(cavpts))
    cavpts  = find(atgetcells(RING, 'FamName', 'CAVITY'));
end
if (isempty(cavpts))
    cavpts  = find(atgetcells(RING, 'FamName', 'cav'));
end
RING = atSetRingProperties(RING,'cavpts',cavpts,'rf_voltage',V0);

LMA=calcLMA(RING,MAoptions);

map_l = LMA.outputs.map_l;
map_h = LMA.outputs.map_h;

%% Touschek lifetime calculation at different RF voltage using previously generated LMA 
Vrf = linspace(1.1*U0, 6*U0, 30); % voltage range
maxLifetime = 0;
Voptimum = 0;

for jj=1:length(Vrf)
    V=Vrf(jj);
    bh=bheight(V,U0,alphac,energy,harm);
    RING = atSetRingProperties(RING,'cavpts',cavpts,'rf_voltage',Vrf(jj));

% Update map_h and map_l corresponding to bucket height 
    map_h_new = min(map_h, bh*ones(length(map_h),1));
    map_l_new = max(map_l, -bh*ones(length(map_l),1));

    LMA_new = LMA;
    LMA_new.outputs.map_h = map_h_new;
    LMA_new.outputs.map_l = map_l_new;

    LT=calcTLT(RING,TLoptions,MAoptions,'LMA',LMA_new);

    Tlife=LT.outputs.TL/3600;

    if Tlife > maxLifetime
        maxLifetime = Tlife;
        Voptimum = V;
    else
        % If Tlife decreases after a certain point, break the loop
        if Tlife < maxLifetime
            break;
        end
    end
end

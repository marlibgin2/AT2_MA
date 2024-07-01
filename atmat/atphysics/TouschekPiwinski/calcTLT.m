function TL = calcTLT(varargin)
% Wrapper for Touschek lifetime calculation. Calculates, if necessary, 
% linear optics parameters; expects a full ring as input, but may choose 
% only one achromat for the Local Momentum Aperture (LMA) calculations
%
%% Inputs
% Mandatory arguments
% RING : AT2.0 lattice cell array.
%
% ***************************************************
% TLoptions: structure with Touschek lifetime calculation options,
%            containing the following fields (if empty, default values are set) : 
%
% TLoptions.Ib : current per bunch [A],. default = 0.5/176
% TLoptions.integrationmethod : see calcTLT_raw, default =
%                                               'integral'
% TLoptions.abstol: absolute tolerance, default  = 1.0e-16
% TLoptions.reltol: relative tolerance, default  = 1.0e-16
% TLoptions.nperiods: number of periods contained in the input lattice, default = 20
% TLoptions.LMAperiods: number of periods for which to calculate the LMA, default = 1;
% TLoptions.kcoupl: coupling ratio, default = 0.025
%
% ************************************
% MAoptions: structure with local momentum aperture calculation options,
%            containing the following fields (if empty, default values are set)
%
% MAoptions: structure with local momentm aperture calculation
%                         options, with fields:
%
% MAoptions.lmafams: cell array of strings with names of magnet families at which LMA
%                     is to be calculated. If = 'all' then all non-zero length elements are included           = 'all';
% MAoptions.stepfam: specifies only one every stepfam elements are included
% MAoptions.deltalimit: maximum momentum deviation to be searched. Used to establish the rf bucket height.
% MAoptions.initcoord: initial coordinates [x0 y0]
% MAoptions.delta: initial guess for momentum aperture 
% MAoptions.deltastepsize: step size for LMA search
% MAoptions.splits: number of iterations of step division
% MAoptions.split_step_divisor: factor to reduce step size at each iteration
% MAoptions.nturns: number of turns. If nan then number of turns
%                             is chosen as 1.2/Qs, default = 1024
% MAoptions.S0max: maximum longitudinal position at which to calculate LMA [m], 
%                            default if RING is available = findspos(cLOptions.RING,length(cLOptions.RING)+1)/20;
%                            default if RING is not available = 528/20
% MAoptions.S0min: minimum longitudinal position at which to
%                            calculate LMA [m], default = 0.0
% Optional arguments
%   verbose     :defines level of verbose output, default=0, i.e. no output 
%
%% Outputs
% TL : structure with the fields
%   TL.inputs echoes the input
%   TL.inputs.RING 
%
%   TL.inputs.TLoptions
%   TL.inputs.MAoptions
%
%   TL.outputs structure with fields

%   TL.outputs.LMA: local momentum aperture structure  with fields
%
%   TL.outputs.LMA.outputs lists calculation results with the following fields
%   TL.outputs.LMA.outputs.Spos  : (1Xnspos) array of longitudinal positions where LMA is calcualated [m]
%   TL.outputs.LMA.outputs.map_l : (1XnSpos) array of negative LMAs
%   TL.outputs.LMA.outputs.map_h : (1XnsPos) array of positive LMAs
%   TL.outputs.LMA.outputs.PeriodDev : if periodicity check was requested, this contains the
%                           deviation from periodicity
%   TL.outputs.LMA.outputs.MAoptions: LMA calculation options
%   TL.outputs.LMA.outputs.telapsed: elapsed calculation time [s]

%   TL.outputs.TL Touschek lifetime [s]
%   TL.outputs.telapsed: elapsed calculation time [s]
%
%% Usage Examples
%
% LT = calcTLT(RING,TLoptions,MAoptions,'verbose',2);

%% History
% PFT 2024/03 first version
% PFT 2024/06/23 : restructured input/output and improved documentation
%
%% Input argument parsing
[RING,TLoptions,MAoptions]= getargs(varargin,[],[],[]);

if isempty(TLoptions)
    TLoptions.Ib = 0.5/176;
    TLoptions.integrationmethod = 'integral';
    TLoptions.abstol = 1.0e-16;
    TLoptions.reltol = 1.0e-16;
    TLoptions.nperiods = 20;
    TLoptions.LMAperiods = 1;
    TLoptions.kcoupl = 0.025;
end

verboselevel  = getoption(varargin,'verbose',0);
Ib = getoption(varargin,'Ib',TLoptions.Ib);
integrationmethod = getoption(varargin,'integrationmethod',TLoptions.integrationmethod);
abstol = getoption(varargin, 'AbsTol', TLoptions.abstol);
reltol = getoption(varargin, 'Relol', TLoptions.reltol);
nperiods = getoption(varargin,'nperiods',TLoptions.nperiods);
LMAperiods = getoption(varargin,'LMAperiods',TLoptions.LMAperiods);
kcoupl = getoption(varargin,'kcoupl',TLoptions.kcoupl);
   
TLoptions.Ib=Ib;
TLoptions.integrationmethod = integrationmethod;
TLoptions.abstol = abstol;
TLoptions.reltol = reltol;
TLoptions.nperiods = nperiods;
TLoptions.LMAperiods = LMAperiods;
TLoptions.kcoupl = kcoupl;

if (isempty(MAoptions))
   MAoptions.lmafams='all';
   MAoptions.stepfam=1;
   MAoptions.deltalimit=0.1;
   MAoptions.initcoord=[0.0 0.0 0.0 0.0 0.0 0.0]';
   MAoptions.delta=0.01;
   MAoptions.deltastepsize=0.001;
   MAoptions.splits=10;
   MAoptions.split_step_divisor=2;
   MAoptions.nturns=500;
   MAoptions.S0max=528/20;
   MAoptions.S0min=0.0;
end

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

%% Calculates Touschek lifetime

tstart=tic;
if (verboselevel>0)
    fprintf('**** \n');
    fprintf('%s calcTLT: starting - calculating atsummary \n', datetime);
end
energy=RING{1}.Energy;
Periodicity=RING{1}.Periodicity;
circumference = findspos(RING,length(RING)+1)*Periodicity;
period = circumference/nperiods;
ats = atsummary(RING);
emitx = ats.naturalEmittance/(1+kcoupl);
emity = ats.naturalEmittance*kcoupl/(1+kcoupl);
sigp  = ats.naturalEnergySpread;
sigs  = ats.bunchlength;
deltalimit = ats.energyacceptance;

S0min = 0.0;
S0max = LMAperiods*period;
if (verboselevel>0)
    fprintf('%s calcTLT: calculating LMA \n', datetime);
end
LMA = calcLMA(RING,MAoptions,'S0max',S0max,'S0min',S0min,...
                  'deltalimit',deltalimit,'verbose', verboselevel-1);

map_l = LMA.outputs.map_l;
map_h = LMA.outputs.map_h;
Ipos  = LMA.outputs.Ipos;

lmap=cat(2,map_h,map_l);
[~,lindata] = atlinopt4(RING,Ipos);
for i=1:length(lindata)
    lindata(i).Length = RING{Ipos(i)}.Length;
end

if (verboselevel>0)
    fprintf('%s calcTLT: calculating TLT \n', datetime);
end

TousLT = calcTLT_raw(lindata,lmap,'Ib',Ib,'circumference',circumference,'energy',energy,'emitx',emitx,'emity',emity,...
          'sig',sigp,'sigs',sigs,'abstol',abstol,'reltol', reltol, 'integrationmethod', integrationmethod,'verbose',verboselevel-1);

%% Collects output structure data
TL.inputs.RING =RING;

TL.inputs.MAoptions=MAoptions;
TL.inputs.TLoptions=TLoptions;
TL.outputs.LMA = LMA;
TL.outputs.TL = TousLT;
TL.outputs.telapsed=toc(tstart);
if (verboselevel>0)
    fprintf('%s calcTLT: TLT calculation complete \n', datetime);
end


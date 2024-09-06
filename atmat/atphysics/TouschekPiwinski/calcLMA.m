function LMA = calcLMA(varargin)
% Calculates and plots Local Momentum Acceptance. Tracking can be 6d or 4d
% as defined by the input lattice. This is intended as a higher level 
% wrapper function that in turn calls the lower level functions
% "calcLMA_raw" and "plotLMA"
%
%% Inputs 
% Mandatory input arguments
% RING : AT2 lattice array
% MAoptions :Structure containing the following fields:
%  MAoptions.lmafams: cell array of strings with names of magnet families at which LMA
%                     is to be calculated. If = 'all' then all non-zero length elements are included
%  MAoptions.stepfam: specifies only one every stepfam elements are included
%  MAoptions.deltalimit: maximum momentum deviation to be searched. Used to establish the rf bucket height.
%  MAoptions.initcoord: initial coordinates [x0 x0p y0 x0p delta z0]'
%  MAoptions.delta: initial guess for momentum aperture 
%  MAoptions.deltastepsize: step size for LMA search;
%  MAoptions.splits : number of iterations of step division
%  MAoptions.split_step_divisor: factor to reduce step size at each iteration
%  MAoptions.nturns: number of turns. if nan use nsynchT synchrotron periods
%  MAoptions.nsyncT: number of synchrotron periods to be used if nturns is nan
%  MAoptions.S0max: maximum longitudinal position at which to calculate LMA
%  MAoptions.S0min: minimum longitudinal position at which to calculate LMA
%  
%            If MAoptions = [], hard-coded defaults are used.
%            Values of MAoptions fields are overridden if given explicitly 
%            as input in the form ('parameter', value)
% 
% Optional input parameters
%
% all fields in MAoptions
%
% nperiods: number of periods contained in the input lattice. Used to 
%           determine the period length for periodicity checks. 
%           RING is assumed to contain the whole ring in this case, 
%           default = 20
%
% verbose : defines level of verbose output, default=0, i.e. no output
%
% Optional flags
% plot : plots LMA
% checkperiodicity: calculates deviation from perididicity i.e,
%                   max(abs(LMA(s+L*i)-LMA(s)) where L is the period length
%                   and i varies from 1 to the number of periods
%
%% Outputs 
% Structure with the following fields
% LMA.inputs echoes the input parameters
%   LMA.inputs.RING : inut lattice
%
% LMA.outputs lists calculation results with the following fields
%   LMA.outputs.Spos  : (1Xnspos) array of longitudinal positions where LMA is calcualated [m]
%   LMA.outputs.map_l : (1XnSpos) array of negative LMAs
%   LMA.outputs.map_h : (1XnsPos) array of positive LMAs
%   LMA.outputs.PeriodDev : if periodicity check was requested, this contains the
%                           deviation from periodicity
%   LMA.outputs.MAoptions: LMA calculation options
%   LMA.outputs.telapsed: elapsed calculation time [s]
%
%% Usage examples
% LMA = calcLMA(RING_a1,MAoptions,'nperiods',20,'S0min',0.0,'S0max',528,'checkperiodicity');
%
% LMA = calcLMA(RING,MAoptions,'nperiods',20,'lmafams','all','verbose',1)
%
% LMA = calcLMA(RING,[],'lmafams',{'S5_a1'},'nturns',nan);
%

%% History
% PFT 2024/03/08. 
% PFT 2024/06/16: changed output into a structure
% PFT 2024/07/12: corrected naming of PeriodDev field in output Structure
% PFT 2024/07/30: changed initcoord default for handling of nan as initial
%                 phase
% PFT 2024/09/06: added handling of nsyncT in MAoptions.

%% Input argument parsing
[RING,MAoptions] = getargs(varargin,[],[]);
if (isempty(MAoptions))
    MAoptions.nperiods=20;
    MAoptions.lmafams='all';
    MAoptions.stepfam=1;
    MAoptions.stepfam=1;
    MAoptions.deltalimit=0.3;
    MAoptions.initcoord=[0 0 0 0 0 nan]';
    MAoptions.delta=0.01;
    MAoptions.deltastepsize=0.1;
    MAoptions.splits=10;
    MAoptions.split_step_divisor=2;
    MAoptions.nturns=nan;
    MAoptions.nsyncT=3;
    MAoptions.S0max=528/20;
    MAoptions.S0min=0.0;
end
plotf              = any(strcmpi(varargin,'plot'));
checkperiodicityf  = any(strcmpi(varargin,'checkperiodicity'));
verboselevel       = getoption(varargin,'verbose',0);
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
if isfield(MAoptions,'nsyncT')
    nsyncT       = getoption(varargin,'nsyncT',MAoptions.nsyncT);
else
    nsyncT       = 3;
end
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
MAoptions.nyncT              = nsyncT;
MAoptions.S0max              = S0max;
MAoptions.S0min              = S0min;

%% Locate points at which LMA is to be calculated
tstart = tic;
if (strcmpi(lmafams,'all'))
    Ipos = (1:size(RING,1))'; %(first element is assumed to be RingPAram)
else
    Ipos = find(atgetcells(RING,'FamName', lmafams{1}));
    for i=2:size(lmafams,1)
        Ipos = cat(1,Ipos,find(atgetcells(RING,'FamName',lmafams{i})));
    end
    Ipos = sort(Ipos);
end
Ipos = Ipos(1:stepfam:size(Ipos,1));
Spos = findspos(RING,Ipos);

if(not(isnan(S0max))&&not(isnan(S0min)))
    Ipos = Ipos(intersect(find(Spos<=S0max)',find(Spos>=S0min)'));
end
Ipos = Ipos(find(atgetfieldvalues(RING, Ipos, 'Length'))); % makes sure only non-zero length elements are included
Spos = findspos(RING,Ipos);
nSpos = numel(Spos);

%% Calculate LMA
if (verboselevel>0)
    fprintf('**** \n');
    fprintf('%s calcLMA: LMA at %3d points \n', datetime, length(Spos));
end

if (isnan(nturns))
    if (verboselevel>0)
        fprintf('%s calcLMA: nturns=nan, calculating atsummary \n', datetime);
    end
    ats=atsummary(RING);

    nturns = round(nsyncT/ats.synctune);
    MAoptions.nturns=nturns;
    if (verboselevel>0)
        fprintf('%s calcLMA: nturns = %3d \n', datetime, nturns);
    end
end

if (verboselevel>0)
        fprintf('%s calcLMA: calculating LMA \n', datetime);
end

if (not(isnan(nturns)))
    [map_l,map_h]=calcLMA_raw(RING,Ipos,...
               'deltalimit',deltalimit, ...
               'initcoord', initcoord,...
               'delta', delta,...
               'deltastepsize',deltastepsize,...
               'splits',splits,...
               'split_step_divisor',split_step_divisor,...
               'nturns',nturns,...
               'verbose',verboselevel-2);
else
    fprintf('%s calcLMA: nturns not defined \n', datetime);
    map_l=zeros(nSpos,1);
    map_h=zeros(nSpos,1);
end
telapsed=toc(tstart);

%% Checks periodicity
if(checkperiodicityf) % assumes a full ring is the input and that S0max (if given) is larger or equal to the circumference
    if (RING{1}.Periodicity>1)
        fprintf('%s CalcPlotMA Warning: Not enough data for periodicity check \n', datetime);
        return
    end
    Circ    = findspos(RING,length(RING)+1);
    Period  = Circ/nperiods;
    Iperiod = find(Spos<=Period)';
        
    nelemperperiod = size(Iperiod,1);
    map_h_max = zeros(nelemperperiod,1);
    map_h_min = zeros(nelemperperiod,1);
    map_l_max = zeros(nelemperperiod,1);
    map_l_min = zeros(nelemperperiod,1);
    for i=1:nelemperperiod
        S=Spos(Iperiod(i,1));
        for j=2:nperiods
            k = find(abs(Spos-(S+(j-1)*Period))<=1E-4,1);
            if (not(isempty(k)))
                Iperiod(i,j) = k(1,1);
            else
                fprintf('%s CalcPlotMA:Warning - could not find corresponding point at S = %10.5f in achromat %3d %3d\n',datetime, S, j, i);
            end
        end
    
        map_h_max(i)=max(map_h(Iperiod(i,1:nperiods)));
        map_h_min(i)=min(map_h(Iperiod(i,1:nperiods)));
        map_l_max(i)=max(map_l(Iperiod(i,1:nperiods)));
        map_l_min(i)=min(map_l(Iperiod(i,1:nperiods)));
    end
    diff_h = map_h_max-map_h_min;
    diff_l = map_l_max-map_l_min;
    Periodicity_h = max(diff_h);
    Periodicity_l = max(diff_l);
    PeriodDev=[Periodicity_l, Periodicity_h]*100;
    fprintf('%s Periodicity Deviation low  = %8.4f %% \n', datetime, Periodicity_l*100);
    fprintf('%s Periodicity Deviation high = %8.4f %% \n', datetime, Periodicity_h*100);
else
    PeriodDev=NaN;
end
%% Collects output structure data
LMA.inputs.RING=RING;
LMA.outputs.MAoptions=MAoptions;
LMA.outputs.Spos  = Spos;
LMA.outputs.map_l = map_l;
LMA.outputs.map_h = map_h;
LMA.outputs.Ipos  = Ipos;
LMA.outputs.PeriodDev = PeriodDev;
LMA.outputs.telapsed = telapsed;

%% Plots LMA
if (plotf)
    plotLMA(LMA)
end


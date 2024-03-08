function [map_l,map_h, Spos, Ipos, PeriodDev] = CalcPlotMA(varargin)
% Calculates and Plots Local Momentum Acceptance
% 
%% Usage examples
% [map_l,map_h, Spos, Ipos, PeriodDev] = CalcPlotMA(RING_a1,LatticeOptData,'nperiods',20,'S0min',0.0,'S0max',528,'checkperiodicity');
%
% [map_l,map_h, Spos, Ipos, ~] = CalcPlotMA(RING,LatticeOptData,'nperiods',20,'lmafams','all')
%
% [map_l,map_h, Spos, Ipos, ~] = CalcPlotMA(RING,[],'lmafams','all','nturns',nan);
%
%% Mandatory input arguments
% RING : AT2 lattice array
% LatticeOptData : Structure containing various default values. If input = [], hard-coded defaults are used.
% 
%% Optional input parameters
%
% nperiods: number of periods - used to determine the period length for
%           periodicity checks. RING is assumed to contain the whole ring in this case
% lmafams: cell array of strings with names of magnet families at which LMA
%          is to be calculated. If = 'all' then all non-zero lenmgth
%          elements are included
% stepfam: specifies only one every stepfam elements are included
% deltalimit: maximum momentum deviation to be searched. Used to establish the rf bucket height.
% initcoord: initial coordinates [x0 x0p y0 x0p delta z0]'
% delta: initial guess for momentum aperture 
% deltastepsize: step size for LMA search;
% splits : number of iterations of step division
% split_step_divisor: factor to reduce step size at each iteration
% nturns: numbr of turns. If nan then number of turns is chosen as 1.2/Qs
%                         this is handled by the momentum:aperture_at
%                         function
% S0max: maximum longitudinal position at which to calculate LMA
% S0min: minimum longitudinal position at which to calculate LMA
%
%% Option flags
% plot : plots LMA
% checkperiodicity: calcualtes deviation from perididicityt
% verbose: produces verbose output
%
%% Output parameters
% map_l : negative LMA
% map_h : positive LMA
% Spos  : longitudinal positions where LMA is calcualated [m]
% PeriodDev : if periodicity check was requested, this contains the
%            deviation from periodicity

%% Input argument parsing
[RING,LatticeOptData] = getargs(varargin,[],[]);
if (isempty(LatticeOptData))
    LatticeOptData.sext_fams='all';
    LatticeOptData.stepfam_lma=1;
    LatticeOptData.stepfam_lma=1;
    LatticeOptData.deltalimit_lma=0.01;
    LatticeOptData.initcoord_lma=[0.0 0.0 0.0 0.0 0.0 0.0]';
    LatticeOptData.delta_lma=0.01;
    LatticeOptData.deltastepsize_lma=0.001;
    LatticeOptData.splits_lma=10;
    LatticeOptData.split_step_divisor_lma=2;
    LatticeOptData.nturns_lma=500;
    LatticeOptData.S0max_lma=528/20;
    LatticeOptData.S0min_lma=0.0;
end
plotf              = any(strcmpi(varargin,'plot'));
checkperiodicityf  = any(strcmpi(varargin,'checkperiodicity'));
verbosef           = any(strcmpi(varargin,'verbose'));
nperiods           = getoption(varargin,'nperiods',20);
lmafams            = getoption(varargin,'lmafams',LatticeOptData.sext_fams);
stepfam            = getoption(varargin,'stepfam',LatticeOptData.stepfam_lma);
deltalimit         = getoption(varargin,'deltalimit',LatticeOptData.deltalimit_lma);
initcoord          = getoption(varargin,'initcoord',LatticeOptData.initcoord_lma);
delta              = getoption(varargin,'delta',LatticeOptData.delta_lma);
deltastepsize      = getoption(varargin,'deltastepsize',LatticeOptData.deltastepsize_lma);
splits             = getoption(varargin,'splits',LatticeOptData.splits_lma);
split_step_divisor = getoption(varargin,'split_step_divisor',LatticeOptData.split_step_divisor_lma);
nturns             = getoption(varargin,'nturns',LatticeOptData.nturns_lma);
S0max              = getoption(varargin,'S0max', LatticeOptData.S0max_lma);
S0min              = getoption(varargin,'S0min', LatticeOptData.S0min_lma);


%% Locate points at which LMA is to be calculated
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

%% Calculate LMA
if (verbosef)
    fprintf('**** \n');
    fprintf('%s CalcPlotMA : Calculating LMA at %3d points \n', datetime, length(Spos));
    tic;
end

if (isnan(nturns))
    ats=atsummary(RING);
    nturns = 1.2/ats.synctune;
end
[map_l,map_h]=MomAperture_allRing_par(RING,Ipos,...
               'deltalimit',deltalimit, ...
               'initcoord', initcoord,...
               'delta', delta,...
               'deltastepsize',deltastepsize,...
               'splits',splits,...
               'split_step_divisor',split_step_divisor,...
               'nturns',nturns);
if (verbosef)
    toc;
end

%% Plots LMA
if (plotf)
    figure;plot(Spos, map_l*100, '-o');hold;plot(Spos,map_h*100,'o-');
    xlabel('S[m]');
    ylabel('Local Momentum Aperture [%]');
    grid on;
    ylim([-15,15]);
end

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
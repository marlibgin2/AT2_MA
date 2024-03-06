function [map_l,map_h, Spos, Ipos, PeriodDev] = CalcPlotMA(varargin)
% Calculates and Plot MA
% 
%% Input argument parsing
[RING,LatticeOptData] = getargs(varargin,[],[]);
plotf              = any(strcmpi(varargin,'plot'));
checkperiodicityf  = any(strcmpi(varargin,'checkperiodicity'));
verbosef           = any(strcmpi(varargin,'verbose'));
nperiods           = getoption(varargin,'nperiods',20);
lmafams            = getoption(varargin,'lmafams',LatticeOptData.sext_fams);
stepfam            = getoption(varargin,'stepfam',LatticeOptData.stepfam_lma);
deltalimit         = getoption(varargin,'deltalimit',LatticeOptData.deltalimit_lma);
initcoord          = getoption(varargin,'initcoord',LatticeOptData.initcoord_lma);
delta              = getoption(varargin,'delta',LatticeOptData.delta_lma);
precdelta          = getoption(varargin,'precdelta',LatticeOptData.precdelta_lma);
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
    Ipos = find(atgetcells(RING, 'FamName', lmafams{1}));
    for i=2:size(lmafams,1)
        Ipos = cat(1,Ipos,find(atgetcells(RING, 'FamName', lmafams{i})));
    end
    Ipos = sort(Ipos);
end
Ipos = Ipos(1:stepfam:size(Ipos,1));
Spos = findspos(RING,Ipos);

if(not(isnan(S0max))&&not(isnan(S0min)))
    Ipos = intersect(find(Spos<=S0max)',find(Spos>=S0min)');
end
Ipos = find(atgetfieldvalues(RING, Ipos, 'Length')); % makes sure only non-zero length elements are included
Spos = findspos(RING,Ipos);

%% Calculate LMA
if (verbosef)
    fprintf('**** \n');
    fprintf('%s CalcPlotMA : Calculating LMA at %3d points \n', datetime, length(Spos));
    tic;
end

[map_l,map_h]=MomAperture_allRing_par(RING,Ipos,...
               'deltalimit',deltalimit, ...
               'initcoord', initcoord,...
               'delta', delta,...
               'precdelta',precdelta,...
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
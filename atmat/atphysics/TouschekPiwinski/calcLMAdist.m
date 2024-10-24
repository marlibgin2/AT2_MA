function LMAdist = calcLMAdist(varargin)
% Calculates and plots the Local Momentum Aperture of a number of 
% lattice variants that differ only through the application of a given error model
% (and possibly corresponding corrections)
% Tracking can be 6d or 4d as defined by the input lattice. 
% This is a higher level wrapper function
% that in turn calls the lower level function "calcLMA_raw"
% 
%% Inputs
% Mandatory arguments
% RING : AT2 lattice array
% ErrorModel: structure generated by the errormodel function : 
%              if empty default is 
%              errormodel_DDRchallenging('gdran',1.0,'magalran',1.0,...
%                                'mulsys',1.0, 'mulran',1.0,...
%                                 'bpmran', 1.0, 'strran', 1.0);
% MAoptions :Structure containing the following fields:
%
%            lmafams: cell array of strings with names of magnet families at which LMA
%                     is to be calculated. If = 'all' then all non-zero length elements are included
%            stepfam: specifies only one every stepfam elements are included
%            deltalimit: maximum momentum deviation to be searched. Used to establish the rf bucket height.
%            initcoord: initial coordinates [x0 x0p y0 x0p delta z0]'
%            delta: initial guess for momentum aperture 
%            deltastepsize: step size for LMA search;
%            splits : number of iterations of step division
%            split_step_divisor: factor to reduce step size at each iteration
%            nturns: number of turns. if nan use nsynchT synchrotron periods
%            nsyncT: number of synchrotron periods to be used if nturns is nan
%            S0max: maximum longitudinal position at which to calculate LMA
%            S0min: minimum longitudinal position at which to calculate LMA
%  
%            If MAoptions = [], hard-coded defaults are used.
%            Values of MAoptions fields are overridden if given explicitly 
%            as input in the form ('parameter', value)
% 
% Optional arguments
% all fields in MAoptions
%
% ERlat: structure generated by the "generate_errlatt" function. If not
%   available or empty this structure is generated in this function by a
%   call to "generate_errlatt".
% corrorb: if true, perform orbit correction
% corrtun: if true, perform tune correction
%
% verbose : defines level of verbose output, default=0, i.e. no output
%
% Optional flags
% plot : plots LMAdist
% fulloutput : includes detailed data all data on each one of seeds.


%% Outputs
% Structure with the following fields
% LMAdist.inputs echoes the input parameters
%   LMAdist.inputs.RING : input lattice
%   LMAdist.inputs.ErrorModel 
%   LMAdist.inputs.nseeds 
%   LMAdist.inputs.tunfams 
%   LMAdist.inputs.nittune 
%   LMAdist.inputs.TolTune 
%   LMAdist.inputs.frac 
%   LMAdist.inputs.corrorbf : correct orbit flag
%   LMAdist.inputs.corrtunf : correct tune flag
%
% LMAdist.outputs lists calculation results with the following fields
%   LMAdist.outputs.desc : : datetime + input description
%   LMAdist.outputs.MAoptions: LMA calculation options
%   LMAdist.outputs.Spos  : (1XnSpos) array of longitudinal positions where LMA is calculated [m]
%   LMAdist.outputs.Ipos  : (1XnSpos) array of indices to the elements at
%                           which LMA is to be calculated
%   LMAdist.outputs.map_l : (nseeds+1xnSpos) array of negative LMA for all seeds
%   LMAdist.outputs.map_h : (nseeds+1xnSpos) array of positive LMA for all seeds
%   LMAdist.outputs.map_l_av: (1Xnspos) array of average negative LMA
%   LMAdist.outputs.map_h_av: (1Xnspos) array of average positive LMA
%   LMAdist.outputs.map_l_std: (1Xnspos) array of standard deviations of  negative LMA
%   LMAdist.outputs.map_h_std: (1Xnspos) array of standard deviations of  positive LMA
%
%   LMAdist.outputs.orb0_stds: (6Xnseeds+1) array of close orbit standadr
%                           deviations for perturbed lattices before 
%                           correction. First point is the unperturbed
%                           lattice.
%   LMAdist.outputsorb_stds  = (6Xnseeds+1) array of close orbit standard
%                           deviations for perturbed lattices after 
%                           correction. First point is the unperturbed
%                           lattice.
% Note: the outptus (RINGe,raparae,Itunese and Ftunese) below are empty 
%       unless the 'fulloutput' option is on
%   LMAdist.outputs.RINGe:  (nseeds+1Xsize of RING) cell array of perturbed 
%                        lattices after correction. first is the
%                        unperturbed lattice.
%   LMAdist.outputs.rparae: (nseeds+1X1) cell array of atsummaries for
%                        perturbed lattices after correction. 
%                        First is the unperturbed lattice.
%   LMAdist.outputs.Itunese:(nseeds+1X1) cell array of tunes for the 
%                         perturbed lattices before correction
%                         First is the unperturbed lattice.
%   LMAdist.outputs.Ftunese: (nseeds+1X1) cell array of tunes for the 
%                         perturbed lattices after correction
%                         First is the unperturbed lattice.
%   LMAdist.outputs.telapsed : elapsed calculation time [s]
%
%% Usage examples
% LMAdist = calcLMAdist(RING,ErrorModel,MAoptions,'plot','corrorb','verbose', 1);
% DAdist = calcDAdist(RING,ErrorModel,[],'nturns',1024,'corrorb');
% calcDAdist(RING,ErrorModel,[],'nturns',1024,'nseeds',10,'plot','corrorb');
% DAdist = calcDAdist(RING,ErrorModel,[],'nturns',810,'corrorb','corrtun','tunfams',{'Q1','Q2'},'frac',0.5);
% DAdist = calcDAdist(RING,ErrorModel,[],'mode,'xydp','nturns',810,'corrorb','corrtun','tunfams',{'Q1','Q2'});

%% History
% PFT 2024/06/16: first version
% PFT 2024/06/26: restructuring and documentation
% PFT 2024/07/23: better handling and flagging of unstable lattices
% SJ  2024/07/29: replaced code for error generation and correction by
%                 call to "generate_errlatt".
%                 added possibility of input of ERlat structure
%                 added possibility of fixing the response matrix for all
%                 seeds.
% PFT 2024/07/30: changed initcoord default for handling of nan as initial
%                 phase
% PFT 2024/08/05: bug fix - incorrect initialization of output vectors
%                 if the number of seeds was larger than the default (10)
% PFT 2024/08/06: bug fix - incorrect extraction of rms orbits from ERlat.
% PFT 2024/08/07: update of MAoptions.nturns when input nturns is nan.
% PFT 2024/09/06: added handling of nsyncT in MAoptions.
%
%% Input argument parsing
[RING,ErrorModel,MAoptions] = getargs(varargin,[],[],[]);
if (isempty(ErrorModel))
    ErrorModel=errormodel_DDRchallenging('gdran',1.0,'magalran',1.0,...
                                         'mulsys',1.0, 'mulran',1.0,...
                                         'bpmran', 1.0, 'strran', 1.0);
end

ERlat =  getoption(varargin,'ERlat',struct());

if (isempty(MAoptions))
    MAoptions.lmafams='all';
    MAoptions.stepfam=1;
    MAoptions.stepfam=1;
    MAoptions.deltalimit=0.3;
    MAoptions.initcoord=[0 0 0 0 0.0 nan]';
    MAoptions.delta=0.01;
    MAoptions.deltastepsize=0.001;
    MAoptions.splits=10;
    MAoptions.split_step_divisor=2;
    MAoptions.nturns=nan;
    MAoptions.nsyncT=3;
    MAoptions.S0max=528/20;
    MAoptions.S0min=0.0;
end

plotf              = any(strcmpi(varargin,'plot'));
plotorbrmsf        = any(strcmpi(varargin,'plotorbrms'));
corrorbf           = getoption(varargin,'corrorb',true);
corrtunf           = getoption(varargin,'corrtun',true);
fulloutputf        = any(strcmpi(varargin,'fulloutput'));
useORM0f           = getoption(varargin,'useORM0',true);
verboselevel       = getoption(varargin,'verbose',0);
desc               = getoption(varargin,'desc','calcLMAdist:');
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

MAoptions.lmafams    = lmafams;
MAoptions.stepfam    = stepfam;
MAoptions.deltalimit = deltalimit;
MAoptions.initcoord  = initcoord;
MAoptions.delta      = delta;
MAoptions.deltastepsize =  deltastepsize;
MAoptions.splits        = splits;
MAoptions.split_step_divisor = split_step_divisor;
MAoptions.nturns             = nturns;

MAoptions.nsyncT             = nsyncT;
MAoptions.S0max              = S0max;
MAoptions.S0min              = S0min;

nseeds           = getoption(varargin,'nseeds',10);
tunfams          = getoption(varargin,'tunfams',{'Q1_b3','Q2_b3'});
nittune          = getoption(varargin,'nittune',10); % max n. of iterations for tune matching
TolTune          = getoption(varargin,'TolTune',1E-3); % tolerance for tune matching
frac             = getoption(varargin,'frac',1.0); % fraction for quad change in each tune fit iteration

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
nSpos = numel(Spos);

%% Preamble

if (verboselevel>0)
    fprintf('**** \n');
    fprintf('%s CalcLMAdist: Starting LMA distribution calculation at %3d points \n', datetime, length(Spos));
end


if (verboselevel>0)
    fprintf('%s CalcLMAdist: calculating atsummary \n', datetime);
end
try
    rpara=atsummary(RING);
    Itunes = rpara.Itunes;
catch ME
     fprintf('%s calcLMAdist: Error \n', datetime);
     fprintf('Error message was:%s \n',ME.message);
     return
end
   
if (isnan(nturns))

    nturns = round(nsyncT/rpara.synctune);
    MAoptions.nturns=nturns;
    if (verboselevel>0)
        fprintf('%s calcLMAdist: nturns = %3d \n', datetime, nturns);
    end
end

%% Generates or reads lattices with errors
if (verboselevel>0)
    fprintf('*** \n');
    fprintf('%s calcLMAdist: starting perturbed lattice calculations and corrections \n', datetime);
end
tstart= tic;
if (verboselevel>0)
    fprintf('*** \n');
    fprintf('%s calcLMAdist: starting perturbed lattice calculations and corrections \n', datetime);
end

if (isempty(fields(ERlat)))
    ERlat = generate_errlatt(RING,ErrorModel,'tunfams',tunfams, ...
    'nseeds',nseeds,'nittune', nittune, 'TolTune', TolTune,...
    'frac', frac, 'useORM0', useORM0f, 'verbose', verboselevel-1);
else
    if (verboselevel>0)
        fprintf('%s calcLMAdist: using previously calculated corrected lattices \n', datetime);
    end
    nseeds    = ERlat.inputs.nseeds;
    corrorbf  = ERlat.inputs.corrorbf;
    corrtunf  = ERlat.inputs.corrtunf;
    useORM0f  = ERlat.inputs.useORM0f;
    tunfams   = ERlat.inputs.tunfams;
    nittune   = ERlat.inputs.nittune;
    TolTune   = ERlat.inputs.TolTune; 
    frac      = ERlat.inputs.frac;
end

RINGe        = ERlat.outputs.RINGe;
orb0_stds    = ERlat.outputs.orb0_stds;
orb_stds     = ERlat.outputs.orb_stds;
stab         = ERlat.outputs.stab;
survivalrate = ERlat.outputs.survivalrate;

if (not(isempty(ERlat.outputs.rparae{1})))
   rpara = ERlat.outputs.rparae{1};
else
     fprintf('%s Error in calcLMAdist: unperturned ring is unstable \n', datetime);
     LMAdist.inputs.RING=RING;
     LMAdist.inputs.ErrorModel=ErrorModel;
     LMAdist.outputs.orb0_stds=orb0_stds;
     LMAdist.outputs.orb_stds=orb_stds;
     LMAdist.outputs.RINGe=RINGe;
     LMAdist.outputs.rparae=rparae;
     LMAdist.outputs.Itunese=Itunesse;
     telapsed=toc(tstart);
     LMAdist.outputs.telapsed=telapsed;
     return
end


map_l     = zeros(nseeds+1,nSpos);
map_h     = zeros(nseeds+1,nSpos);
map_l_av  = zeros(1,nSpos);
map_h_av  = zeros(1,nSpos);
map_l_std = zeros(1,nSpos);
map_h_std = zeros(1,nSpos);
%

rparae    = cell(nseeds+1,1);
Itunese   = cell(nseeds+1,1);
Ftunese   = cell(nseeds+1,1);

%% Calculate LMAs
if (verboselevel>0)
    fprintf('*** \n');
    fprintf('%s starting LMAdist calculations \n', datetime);
end

for i=1:nseeds+1
 if (verboselevel>1)
    fprintf('%s calcLMAdist: seed n. %3d \n', datetime, i-1);
 end

 if (stab(i))
    RINGes=RINGe{i}; 
    [map_l(i,:),map_h(i,:)]=calcLMA_raw(RINGes,Ipos,...
               'deltalimit',deltalimit, ...
               'initcoord', initcoord,...
               'delta', delta,...
               'deltastepsize',deltastepsize,...
               'splits',splits,...
               'split_step_divisor',split_step_divisor,...
               'nturns',nturns,...
               'verbose',verboselevel-2);
 else
     map_l(i,:)=nan;
     map_h(i,:)=nan;
 end
end
map_l_av(:)=mean(map_l(2:nseeds+1,:),'omitnan');
map_h_av(:)=mean(map_h(2:nseeds+1,:),'omitnan');
map_l_std(:)=std(map_l(2:nseeds+1,:),'omitnan');
map_h_std(:)=std(map_h(2:nseeds+1,:),'omitnan');

telapsed=toc(tstart);
%% Collects output structure data

LMAdist.inputs.RING       = RING;
LMAdist.inputs.ErrorModel = ErrorModel;
LMAdist.inputs.nseeds     = nseeds;
LMAdist.inputs.tunfams    = tunfams;
LMAdist.inputs.nittune    = nittune;
LMAdist.inputs.TolTune    = TolTune;
LMAdist.inputs.frac       = frac;
LMAdist.inputs.corrorbf   = corrorbf;
LMAdist.inputs.corrtunf   = corrtunf;


LMAdist.outputs.desc=strcat(sprintf('%s',datetime),' : ', desc);
LMAdist.outputs.MAoptions=MAoptions;
LMAdist.outputs.Spos=Spos;
LMAdist.outputs.Ipos=Ipos;
LMAdist.outputs.map_l=map_l;
LMAdist.outputs.map_h=map_h;
LMAdist.outputs.map_l_av=map_l_av;
LMAdist.outputs.map_h_av=map_h_av;
LMAdist.outputs.map_l_std=map_l_std;
LMAdist.outputs.map_h_std=map_h_std;

LMAdist.outputs.orb0_stds=orb0_stds;
LMAdist.outputs.orb_stds=orb_stds;
if (fulloutputf)
    LMAdist.outputs.RINGe=RINGe;
    LMAdist.outputs.rparae=rparae;
    LMAdist.outputs.Itunese=Itunese;
    LMAdist.outputs.Ftunese=Ftunese;
else
    LMAdist.outputs.RINGe={};
    LMAdist.outputs.rparae={};
    LMAdist.outputs.Itunese={};
    LMAdist.outputs.Ftunese={};
end
LMAdist.outputs.stab         = stab;
LMAdist.outputs.survivalrate = survivalrate;

LMAdist.outputs.telapsed=telapsed;

if(verboselevel>0)
    fprintf('%s LMAdist calculation complete \n', datetime);
end

%% Plots LMA Distribution and rms orbits
if (plotf)
    if (verboselevel>0)
        fprintf('Plotting LMAdist... \n');
        if (plotorbrmsf)
            plotLMAdist(LMAdist,'verbose',verboselevel-1,'plotorbrms');
        else
            plotLMAdist(LMAdist,'verbose',verboselevel-1);
        end
    else
        if (plotorbrmsf)
            plotDAdist(LMAdist,'plotorbrms');
        else
            plotDAdist(LMAdist);
        end
    end
end

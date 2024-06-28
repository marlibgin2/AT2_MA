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
%              errormodel_DDRchallenging('gdran',0.0,'magalran',1.0,...
%                                'mulsys',0.0, 'mulran',1.0,...
%                                 'bpmran', 0.0, 'strran', 0.0);
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
%            nturns: numbr of turns. If nan then number of turns is chosen as 1.2/Qs           
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
% corrorb: if true, perform orbit correction
% corrtun: if true, perform tune correction
%
% verbose : defines level of verbose output, default=0, i.e. no output
%
% Optional flags
% plot : plots LMAdist


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
% PFT 2024/06/16, first version
% PFT 2024/06/26, restructuring and documentation
%% Input argument parsing
[RING,ErrorModel,MAoptions] = getargs(varargin,[],[],[]);
if (isempty(ErrorModel))
    ErrorModel=errormodel_DDRchallenging('gdran',0.0,'magalran',1.0,...
                                         'mulsys',0.0, 'mulran',1.0,...
                                         'bpmran', 0.0, 'strran', 0.0);
end
if (isempty(MAoptions))
    MAoptions.lmafams='all';
    MAoptions.stepfam=1;
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

plotf              = any(strcmpi(varargin,'plot'));
plotorbrmsf        = any(strcmpi(varargin,'plotorbrms'));
corrorbf           = getoption(varargin,'corrorb',true);
corrtunf           = getoption(varargin,'corrtun',true);
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
tstart = tic;
if (verboselevel>0)
    fprintf('**** \n');
    fprintf('%s CalcLMAdist: Starting LMA distribution calculation at %3d points \n', datetime, length(Spos));
end

map_l     = zeros(nseeds+1,nSpos);
map_h     = zeros(nseeds+1,nSpos);
map_l_av  = zeros(1,nSpos);
map_h_av  = zeros(1,nSpos);
map_l_std = zeros(1,nSpos);
map_h_std = zeros(1,nSpos);
%
orb0_stds = zeros(6,nseeds+1);
orb_stds  = zeros(6,nseeds+1);
RINGe     = cell(nseeds+1,1);
rparae    = cell(nseeds+1,1);
Itunese   = cell(nseeds+1,1);
Ftunese   = cell(nseeds+1,1);
stablat   = ones(nseeds+1,1);

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
    nturns = 1.2/rpara.synctune;
end

%% Calculate and correct lattices
if (verboselevel>0)
    fprintf('*** \n');
    fprintf('%s calcLMAdist: starting perturbed lattice calculations and corrections \n', datetime);
end

parfor i=1:nseeds+1
 if (verboselevel>1)
     fprintf('%s calcLMAdist: seed n. %4d \n', datetime, i-1);
 end

 if (i>1)
     RINGe{i}=applyErrorModel(RING,ErrorModel);
     if (corrorbf) 
         try
            if (verboselevel>1)
                fprintf('%s calcLMAdist: Correcting orbit seed n. %3d \n', datetime, i-1);
            end
            [RINGe{i}, orb0, orb] = calcOrb(RINGe{i},'correct');
            for j=1:6
                orb0_stds(j,i)=std(orb0(j,:));
                orb_stds(j,i)=std(orb(j,:));
            end
         catch ME
             fprintf('%s calcLMAdist: Error in orbit correction \n', datetime);
             fprintf('Error message is %s \n', ME.message);
             stablat(i)=0;
         end
     end
     if (corrtunf)
         try
           rparae{i}=atsummary(RINGe{i});
           Itunese{i}=rparae{i}.Itunes;
           if (not(isnan(Itunese{i}(1)))&&not(isnan(Itunese{i}(2))))
              if (verboselevel>1)
                 fprintf('%s calcLMAdist: Fitting tunes from [ %5.3f %5.3f ] to [ %5.3f %5.3f ] seed n. %3d \n',...
                    datetime, Itunese{i}(1),Itunese{i}(2),Itunes(1),Itunes(2), i-1);
              end
              [RINGe{i}, its, penalty_tune, Ftunese{i}] = ...
                      fittuneRS(RINGe{i}, Itunes, tunfams{1}, tunfams{2},...
                      'maxits', nittune,'Tol', TolTune,...
                      'UseIntegerPart',true,'frac',frac,...
                      'verbose',verboselevel-2);
                if (verboselevel>1)
                    fprintf('%s calcLMAdist: Tune fit complete with penalty = %6.2e after %3d iterations seed n. %3d \n', datetime, penalty_tune, its, i-1);
                end
           else
                 fprintf('%s calcLMAdist: Unstable Lattice Tunes = [ %5.3f %5.3f ]  \n',...
                    datetime, Itunese{i}(1), Itunese{i}(2));
           end
         catch ME
             fprintf('%s calcLMAdist: Error in tune correction \n', datetime);
             fprintf('Error message is %s \n', ME.message);
             stablat(i)=0;
         end
     else
         Itunese{i}=[NaN NaN];
         Ftunese{i}=[NaN NaN];
     end
 else
     RINGe{i}=RING;
     rparae{i}=rpara;
     Itunese{i}=Itunes;
     Ftunese{i}=Itunes;
 end
end

%% Calculate LMAs
if (verboselevel>0)
    fprintf('*** \n');
    fprintf('%s starting LMAdist calculations \n', datetime);
end


for i=1:nseeds+1
 if (verboselevel>1)
    fprintf('%s calcLMAdist: seed n. %3d \n', datetime, i-1);
 end

 if (stablat(i))
    [map_l(i,:),map_h(i,:)]=calcLMA_raw(RINGe{i},Ipos,...
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
map_l_av(:)=mean(map_l(2:nseeds+1,:));
map_h_av(:)=mean(map_h(2:nseeds+1,:));
map_l_std(:)=std(map_l(2:nseeds+1,:));
map_h_std(:)=std(map_h(2:nseeds+1,:));

telapsed=toc(tstart);
%% Collects output structure data

LMAdist.inputs.RING = RING;
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
LMAdist.outputs.RINGe=RINGe;
LMAdist.outputs.rparae=rparae;
LMAdist.outputs.Itunese=Itunese;
LMAdist.outputs.Ftunese=Ftunese;

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

function ERlat= generate_errlatt(varargin)
% Generates lattice with errors and orbit/tune corrections
%% Inputs
% Mandatory arguments
% RING: AT2 lattice array - a full 6D-enabled ring  is expected
% ErrorModel: structure generated by the errormodel function : 
%              if input is an empty array ([]), default is 
%              errormodel_DDRchallenging('gdran',0.0,'magalran',1.0,...
%                                'mulsys',1.0, 'mulran',1.0,...
%                                'bpmran',1.0, 'strran',1.0);
% Optional arguments
%   nseeds  : number of seeds, default = 10
%   tunfams : list of magnet families used for ring tun matching, default = {'Q1_b3','Q2_b3'}
%   nittune : number of iterations for tune matching, default = 10
%   TolTune : tolerance for tune matching, default = 1E-3
%   frac    : fraction for quad change in each tune fit iteration, defaut = 1.0
%
%   corrorb   : if true, perform orbit correction, default is true
%   corrtun   : if true, perform tune correction, default is true
%   useORM0   : if true, sets the orbit correction to use the orbit respose
%               matrix for the unperturbed ring for all iterations, 
%               default=true
%   fulloutput: if true, exports also the rparae,Itunese and Ftunese
%               fields, if false, rpare contains only results for the
%               unperturbed lattice
%               the default=false. 
%   verbose : defines level of verbose output, default=0, i.e. no output
%
% Optional flags

%% Outputs
% ERLat structure with the following fields
%
% ERLat.inputs - echoes input parameters
%   ERlat.inputs.nseeds
%   ERlat.inputs.tunfams
%   ERlat.inputs.nittune
%   ERlat.inputs.Toltune
%   ERlat.inputs.frac
%   ERlat.inputs.corrorb
%   ERlat.inputs.corrtun
%   ERlat.inputs.useORM0f
%   ERlat.inputs.fulloutputf

% ERlat.outputs structure with the following fields (first column 
% corresponds to the lattice without errors)
%   ERlat.outputs.RINGe    :(1xnseeds+1) cell array of lattices with errors and corrections
%   ERlat.outputs.rparae   :(1xnseeds+1) cell array of atsumary results
%   ERlat.outputs.Itunese  :(1xnseeds+1) cell array of betatron tunes just after putting errors in machine
%   ERlat.outputs.Ftunese  :(1xnseeds+1) cell array of betatron tunes after the correction of tune 
%   ERlat.outputs.orb0_stds:(1xnsseds+1) array of rms of the closed orbit along the ring before correction 
%   ERlat.outputs.orb_stds :(1xnsseds+1) array of rms of the closed orbit along the ring after correction 
%   ERlat.outputs.stab     :(1xnseeds+1) array of integers = 1  if lattice
%                           is stable, 0 if it is unstable
%   ERlat.outputs.survivalrate : percentage of stable lattices
%   ERlat.outputs.telapsed     : calculation time in s

%% Usage examples
%  ERlat= generate_errlatt(RING,ErrorModel, 'tunfams',{'Q1';'Q2'}, 'verbose', 1);

%% History
% 2024/07/24: SJ, : first version
% 2024/07/25: PFT : restructured inputs, including Errormodel
%                   restructure outputs to ease integration into the
%                   cLatt package
%                   added option to calculate reposnse matrix for the 
%                   unperturbed lattice only
%                   
%
%% Input argument parsing
[RING,ErrorModel] = getargs(varargin,[],[]);
if (isempty(ErrorModel))
    ErrorModel=errormodel_DDRchallenging('gdran',1.0,'magalran',1.0,...
                                         'mulsys',1.0, 'mulran',1.0,...
                                         'bpmran',1.0, 'strran',1.0);
end
corrorbf         = getoption(varargin,'corrorb',true);
corrtunf         = getoption(varargin,'corrtun',true);
verboselevel     = getoption(varargin,'verbose',0);
fulloutputf      = getoption(varargin,'fulloutput',false);
useORM0f         = getoption(varargin,'useORM0',true);

nseeds           = getoption(varargin,'nseeds',10);
tunfams          = getoption(varargin,'tunfams',{'Q1_b3','Q2_b3'});
nittune          = getoption(varargin,'nittune',10); 
TolTune          = getoption(varargin,'TolTune',1E-3); 
frac             = getoption(varargin,'frac',1.0); 

%% preamble
orb0_stds = zeros(6,nseeds+1);
orb_stds  = zeros(6,nseeds+1);
RINGe     = cell(nseeds+1,1);
rparae    = cell(nseeds+1,1);
Itunese   = cell(nseeds+1,1);
Ftunese   = cell(nseeds+1,1);
stab      = ones(1,nseeds+1); 

ERlat.inputs.nseeds      = nseeds;
ERlat.inputs.tunfams     = tunfams;
ERlat.inputs.nittune     = nittune;
ERlat.inputs.TolTune     = TolTune;
ERlat.inputs.frac        = frac;
ERlat.inputs.corrorbf    = corrorbf;
ERlat.inputs.corrtunf    = corrtunf;
ERlat.inputs.useORM0f    = useORM0f;
ERlat.inputs.fulloutputf = fulloutputf;


%% Apply Errors and corrections
if (verboselevel>0)
   fprintf('%s generate_errlat: Calculating unperturbed lattice parameters \n', datetime);
end
tstart=tic;
try
   rpara   = atsummary(RING);
   Itunes  = rpara.Itunes;
   stab(1) = 1;
catch ME
     fprintf('%s Error in generate_errlatt. atsummary of lattice wihtout errors \n', datetime);
     fprintf('Error message was:%s \n',ME.message);
     ERlat.outputs.RINGe   = RINGe;
     ERlat.outputs.rparae  = rparae;
     ERlat.outputs.Itunese = Itunese;
     ERlat.outputs.Ftunese = Ftunese;
     ERlat.outputs.stab=stab;
     telapsed=toc(tstart);
     ERlat.outputs.telapsed=telapsed;
     return
end

if (useORM0f)
    iBPM = findcells(RING,'FamName','BPM');
    if (isempty(iBPM))
        iBPM=findcells(RING,'FamName','mon');
    end
    indHCor=find(atgetcells(RING,'iscorH','H'));
    if (isempty(indHCor))
        indHCor=findcells(RING,'FamName','ch');
    end
    indVCor=find(atgetcells(RING,'iscorV','V'));
    if (isempty(indVCor))
        indVCor=findcells(RING,'FamName','cv');
    end

    ORM = getlinearrespmat(RING,iBPM,indHCor,indVCor);
else
    ORM = [];
end

parfor i=1:nseeds+1
    if (verboselevel)
        fprintf('%s seed n. %4d \n', datetime, i-1);
    end
    if (i>1)
        RINGe{i}=applyErrorModel(RING,ErrorModel);
        if (corrorbf) 
            try
                if (verboselevel>0)
                    fprintf('%s Correcting orbit seed n. %3d \n', datetime, i-1);
                end
                [RINGe{i}, orb0, orb] = calcOrb(RINGe{i},'correct',...
                    'ORM', ORM, 'verbose',verboselevel-1);
                for j=1:6
                    orb0_stds(j,i)=std(orb0(j,:),'omitnan');
                    orb_stds(j,i)=std(orb(j,:),'omitnan');
                end
            catch ME
                fprintf('%s generate_errlatt: Error in orbit correction for seed n. %3d \n', datetime, i-1);
                fprintf('Error message is %s \n', ME.message);
                stab(i)=0;
            end
        end
        if (corrtunf)
            try
                rparae{i}=atsummary(RINGe{i});
                Itunese{i}=rparae{i}.Itunes;
                if (not(isnan(Itunese{i}(1)))&&not(isnan(Itunese{i}(2))))
                    if (verboselevel>0)
                        fprintf('%s Fitting tunes from [ %5.3f %5.3f ] to [ %5.3f %5.3f ] seed n. %3d \n',...
                        datetime, Itunese{i}(1),Itunese{i}(2),Itunes(1),Itunes(2), i-1);
                    end
                    [RINGe{i}, its, penalty_tune, Ftunese{i}] = ...
                      fittuneRS(RINGe{i}, Itunes,tunfams{1}, tunfams{2},...
                      'maxits', nittune,'Tol', TolTune,...
                      'UseIntegerPart',true,'frac',frac,...
                      'verbose',verboselevel-1);
                    if (verboselevel>0)
                        fprintf('%s Tune fit complete with penalty = %6.2e after %3d iterations seed n. %3d \n', datetime, penalty_tune, its, i-1);
                    end
                else
                    fprintf('%s Unstable Lattice Tunes = [ %5.3f %5.3f ]  \n',...
                        datetime, Itunese{i}(1), Itunese{i}(2));
                    stab(i)=0;
                end
            catch ME
                fprintf('%s generate_errlatt: Error in tune correction for seed n. %3d \n', datetime, i-1);
                fprintf('Error message is %s \n', ME.message);
                stab(i)=0;
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
survivalrate = sum(stab(2:nseeds+1))/nseeds*100;

%% Collects output structure data
ERlat.outputs.RINGe   = RINGe;
ERlat.outputs.rparae{1}  = rparae{1};
if (fulloutputf)
    ERlat.outputs.rparae  = rparae;
    ERlat.outputs.Itunese = Itunese;
    ERlat.outputs.Ftunese = Ftunese;
end
ERlat.outputs.orb0_stds = orb0_stds;
ERlat.outputs.orb_stds  = orb_stds;
ERlat.outputs.stab=stab;
ERlat.outputs.survivalrate = survivalrate;
telapsed=toc(tstart);
ERlat.outputs.telapsed=telapsed;



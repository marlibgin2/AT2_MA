function rp=ExMOGA(varargin)
%% General description 
% Examines result of MOGA Scan:
% Produces optics function plots,
% Matches dispersion, beta functions and achromat tunes.
% Calculates and plots dynamic aperture. 
% Performs a ring tune scan.
% Saves OPA input file
%
%% Inputs
% Mandatory arguments:
%           MOGAResults: structure containing MOGA results (see MOGA.m)
%           index     : index of line in "Values" in MOGAResults structure from which 
%                       data will be retrieved. If negative the line number
%                       is used directly. Convenient to keep track of an
%                       individual regardless of the order of sorting of the
%                       ParetoFront array.
%
% Optional arguments input with syntax ('ParameterName', parameter value)
%           betas0    : 1x2 matrix of desired long straight hor/vert beta
%                       functions.
%           tunesuc   : 1x2 matrix of desired unit cell fractional tunes. (isolated cell)    
%           tunes     : 1x2 matrix of desired full ring (20 achromats) tunes.
%           IQrange   : integer tune ranges [Qxmin Qxmax Qymin Qymax]=[46 47 15 18]
%           qrange    : fractional tune range [qxmin qxmax qymin qymax]=[0.1 0.4 0.1 0.4]
%           Npq          = getoption(varargin,'Npq',[10 10]);
%           LatticeOptData: structures with fields that set various tolerances 
%                           for fits. If not given uses
%                           MOGAResults.LatticeOptData
%           DAoptions : alternative set of parameters for DA calculation.
%                        if not given , takes the values from the 
%                        MOGAResults.LatticeOptData.DAoptions structure
%           split     : number o slices into whihc to split the lattice for
%                      field and geometry calculations
%           ACHRO_ref : reference achromat or unit cell used to calculate
%                       orbit deviations
%           verbose  : controls level of diagnostic printout, 0= no output
%
% Optional flags: presence in the argument list activates the option
%           plot     : plots lattice functions for the extracted lattice
%                       as well for the descendent version with dispersion and tune fits
%           fituc    : fits the unit cell tunes to specified values
%           fitdisp  : matches dispersion to zero in long straight
%           fitbeta0 : fits the horizontal and vertical beta
%                       functions at centre of long stright
%           fittune  : fits lattice to match dispersion and fit period tunes
%                       all other fits.
%           fitchrom : fits chromaticity.(if chrom_fams is not empty)
%           plotda   : plots calculates and plots the on-energy
%                       dynamic  aperture.
%           tunescan  : performs a scan of full ring tunes.
%           plottunescan : plots the DA tune scan
%           saveOPA  : saves a text file with the decision variables set for OPA run. 
%                       the file is based on a template - There are two
%                       templates, one used in "simplified" mode (dipoles hae a single slice)and the other
%                       in "complete" mode (dipoles are sliced as in the
%                       standard lattice.
%
%% Outputs
%          rp: structure containing two substructures (inputs and outputs).
%          rp.inputs documents all inputs to ExMOGA
%
%          rp.outputs
%          rp.outputs.atsummary=structure obtained from atsummary_fast on
%                               the final latice for one achromat
%          rp.outputs.K0 : inital gradient as given from MOGA  
%          rp.outputs.Knew = final gradients as given by linear lattice
%                            fits done by ExMOGA
%          rp.outputs.tunesuc0: initial unit cell tunes
%          rp.outputs.tunesuc1: final unit cell tunes
%          rp.outputs.index=iindex;
%          rp.outputs.DAoptions = DAoptions
%          rp.outputs.Sc1 :strength of first chromaticty correction sextupole family
%          rp.outputs.Sc2 :strength of second chromaticty correction sextupole family
%          rp.outputs.DA  :dynamic aperture (rd calc)
%          rp.outputs.FG  :field and gradients profile
%
%% Usage Examples
% rp=ExMOGA(MOGAResults,12,'plot','fitchrom','plotda','verbose');
% rp=ExMOGA(MOGAResults,12,'plot','fitchrom','fitdisp','plotda','verbose');
% rp=ExMOGA(MOGAResults,12,'plot','fitchrom','fitdisp','fittune','tunes',[55.13658 16.20460]);
% rp=ExMOGA(MOGAResults,12,'plot','fitchrom','fitdisp','verbose','saveOPA');
% rp=ExMOGA(MOGAResults,12,'plot','fitchrom','fitdisp','verbose','DAoptions',DAoptions);

%% History
% PFT july 2023, first version
% PFT 2024/08/28 : updates to handle MGA done on a unit cell, general
%                  structure and documentation improvements.
% PFT 2024/08/29 : added field and gradient profiles calculations and plots
%                  updated hanbling of verbose level
%                  added calculation of orbit deviations
% PFT 2024/08/30 : updated handling of Kall variables

%% Input argument parsing
[MOGAResults,index]=getargs(varargin,[],1);
plotf          = any(strcmpi(varargin,'plot'));
fitucf         = any(strcmpi(varargin,'fituc'));
fitdispf       = any(strcmpi(varargin,'fitdisp')); 
fitbeta0f      = any(strcmpi(varargin,'fitbeta0')); 
fittunef       = any(strcmpi(varargin,'fittune'));
fitchromf      = any(strcmpi(varargin,'fitchro')); 
plotdaf        = any(strcmpi(varargin,'plotda')); 
tunescanf      = any(strcmpi(varargin,'tunescan'));
plottunescanf  = any(strcmpi(varargin,'plottunscan'));
saveOPAf       = any(strcmpi(varargin,'saveOPA'));

tunesuc        = getoption(varargin,'tunesuc',[3/7 1/7]);
betas0         = getoption(varargin,'betas0',[9.0 2.0]);
tunes          = getoption(varargin,'tunes',[46.20 16.28]);
Qrange         = getoption(varargin,'IQrange',[46 47 15 18]);
qrange         = getoption(varargin,'qrange',[0.1 0.4 0.1 0.4]);
Npq            = getoption(varargin,'Npq',[10 10]);
split          = getoption(varargin,'split',1);
LatticeOptData = getoption(varargin,'LatticeOptData',MOGAResults.LatticeOptData);
DAoptions      = getoption(varargin,'DAoptions',LatticeOptData.DAoptions);
chroms0        = getoption(varargin,'chroms0',LatticeOptData.DAoptions.chroms0);
DAoptions.chroms0 = chroms0;
ACHRO_ref      = getoption(varargin,'ACHRO_ref',{});
verboselevel   = getoption(varargin,'verbose',0);

%% Preamble
if ( (verboselevel>0)||tunescanf)
    tic;
end

isdipole    = LatticeOptData.isdipole;
nvars       = LatticeOptData.nvars;
LAT         = LatticeOptData.ACHRO;
UC          = LatticeOptData.UC;
Nper        = LatticeOptData.Nper;
scan_fams   = LatticeOptData.scan_fams;

RunNumber   = MOGAResults.RunNumber;
lb          = MOGAResults.lb;
ub          = MOGAResults.ub;

Sc1 = nan;
Sc2 = nan;
DAS = [];

%% Checks for backward compatibility
%
if (isfield(LatticeOptData,'optMode'))
    optMode  = LatticeOptData.optMode;
else
    optMode = 'SIMP';
end

if (isfield(LatticeOptData,'revBmode'))
    revBmode  = LatticeOptData.revBmode;
else
    revBmode = 'All';
end

if (isfield(LatticeOptData,'chrom_fams'))
    chrom_fams = LatticeOptData.chrom_fams; % list of sextupole families to be used for chromaticity correction
else
    chrom_fams = {'hsd_rbchro';'hsfi_rbchro'};
end

if (isfield(LatticeOptData,'fitvars'))
    nfitvars    = size(LatticeOptData.fitvars,1);
else
    nfitvars    = size(MOGAResults.fitvars,1);
end

if(isfield(LatticeOptData,'lin_fams'))
    lin_fams = LatticeOptData.lin_fams;
else
    lin_fams = {};
end

nlinfams = size(lin_fams,1);

if(isfield(LatticeOptData,'uctune_fams'))
    uctune_fams = LatticeOptData.uctune_fams;
else
    uctune_fams{1} = 'reversebend_sim';
    if (strcmpi(optMode,'SIMP'))
        uctune_fams{2} = 'sshdip';
    else
        uctune_fams{2} = 'dip';
    end
end
    
if(isfield(LatticeOptData,'disp_fams'))
    disp_fams = LatticeOptData.disp_fams;
else
    switch optMode
        case {'SIMP','COMP'}
            disp_fams{1} = 'reversebendmc1_sim';
        case 'CHRO'
            switch revBmode
                case 'All'
                    disp_fams{1} = 'reversebendmc1_chro';
                case 'U3'
                    disp_fams{1} = 'qfm1_rbu3chro';
            end
    end
end

if(isfield(LatticeOptData,'betaxy0_fams'))
    betaxy0_fams = LatticeOptData.betaxy0_fams;
else
    switch optMode
     case {'COMP','SIMP'}
         betaxy0_fams = {'qfend_sim';'qdend_sim'};
     case 'CHRO'
         betaxy0_fams = {'qfend_rbu3chro';'qdend_rbu3chro'};
    end
end

if(isfield(LatticeOptData,'ringtune_fams'))
    ringtune_fams = LatticeOptData.ringtune_fams;
else
    switch optMode
     case {'COMP','SIMP'}
         ringtune_fams = {'qfend_sim';'qdend_sim'};
     case 'CHRO'
         ringtune_fams = {'qfend_rbu3chro';'qdend_rbu3chro'};
    end
end

if(isfield(LatticeOptData,'lattMode'))
    lattMode = LatticeOptData.lattMode;
else
    switch revBmode
        case "All"
            lattMode= 'c1';
        case "U3"
            lattMode= 'b1';
        otherwise
            fprinf ('unknown deprecated lattice mode %s \n', revBmode)
    end
end

if(isfield(LatticeOptData,'nittune'))
    nittune = LatticeOptData.nittune;
else
    nittune = 1000;
end
if(isfield(LatticeOptData,'nitdisp'))
    nitdisp = LatticeOptData.nitdisp;
else
    nitdisp = 1000;
end
if (isfield(LatticeOptData,'TolTune'))
    TolTune = LatticeOptData.TolTune;
else
    TolTune = 1E-7;
end
if (isfield(LatticeOptData,'TolDisp'))
    TolDisp = LatticeOptData.TolDisp;
else
    TolDisp = 1E-8;
end
if (isfield(LatticeOptData,'RINGGRD'))
    RINGGRD=LatticeOptData.RINGGRD;
end
if (isfield(LatticeOptData,'ACHROGRD'))
    ACHROGRD=LatticeOptData.ACHROGRD;
end

%% Parameters for dynamic aperture calculation
% note that these may be different from the parameters at the time
% MOGA was run and recorded in the structure 
% MOGAResults.LatticeOptData.DAoptions
%
% If the argument DAOptions given to the function is an empty matrix, the
% corresponding DAOptions are taken from LatticeOptData structure in 
% MOGAResults. 
%

chroms0      = DAoptions.chroms0; % Target chromaticity for one superperiod or unit cell
TolChrom     = DAoptions.TolChrom;% Chromaticity tolerances
Nitchro      = DAoptions.Nitchro; % Max n. iterations of chromaticity correction


% For backward compatibilty
if (~isfield(DAoptions,'xmaxdas'))
    DAoptions.xmaxdas = 0.007;
end

if (~isfield(DAoptions,'xmindas'))
    DAoptions.xmindas = -0.015;
end

if (~isfield(DAoptions,'ymaxdas'))
    DAoptions.ymaxdas = 0.003;
end
   
if(~(isfield(DAoptions,'XmaxDA')))
    DAoptions.XmaxDA = 0.015;
end

if(~(isfield(DAoptions,'YmaxDA')))
    DAoptions.YmaxDA     = 0.004;
end

XmaxDA = DAoptions.XmaxDA;
YmaxDA = DAoptions.YmaxDA;

if(~(isfield(DAoptions,'DAmode')))
    DAoptions.DAmode = 'border';
end

if (strcmp(DAoptions.DAmode,'grid')&&any(strcmpi(varargin,'DAoptions')))
  
%
% Recalculates X0da and Y0da in case the data in DAoptions 
% does not come from the MOGAResults structure
%
    npdax      = DAoptions.npdax;   % number of grid points in x direction is 2*npdax+1
    npday      = DAoptions.npday;   % number of grid points in y direction is  npday+1
    npDA       = DAoptions.npDA;    % total numbr of grid points
    dx = XmaxDA/npdax; % grid stepsize in x [m]
    dy = YmaxDA/npday;   % grid stepsize in y [m]
    X0da = zeros(npDA,1);  % horizontal coordinates of grid points [m]
    Y0da = zeros(npDA,1);  % vertical coordinates of grid points [m]
    k= 1;
    for i=0:npday 
        for j=0:2*npdax
        X0da(k) = -XmaxDA+dx*j;
        Y0da(k) =  dy*i;
        k=k+1;
        end
    end
    DAoptions.X0da=X0da;
    DAoptions.Y0da=Y0da;
end

%% Set Lattice Parameters
%
if (index<0)
    i=abs(index);
else
    indexcol = size(MOGAResults.ParetoFront,2); 
    %
    % Searches for desired individual in the Pareto Front
    % Last column of ParetoFront table must contain index
    %
    try       
        i=find(MOGAResults.ParetoFront(:,indexcol)==index,1);
        if  isempty(i)
            if (verboselevel>0) 
                fprintf('Index %4d not found in ParetoFront; using i=1', index); 
            end
            i=1;
        end
    catch
        if (verboselevel>0) 
            fprintf('Index column not found in ParetoFront; using i=1');
        end
        i=1;
    end
end

i0   = i;
Ks   = MOGAResults.ParetoFront(i0,1:nvars);
Knew = Ks;
%{
if (isfield(LatticeOptData,'All_fams'))
   All_fams = LatticeOptData.All_fams;
   if (not(strcmpi(optMode,'UC')))
     Idvs = zeros(LatticeOptData.nvars,1);
     for i=1:nvars
       if (isfield(LatticeOptData,'All_famsO'))
          Idvs(i) = find(strcmp(LatticeOptData.All_famsO,scan_fams{i}));
       else
          Idvs(i) = find(strcmp(All_fams,scan_fams{i})); 
       end
     end
     Ilin = zeros(nlinfams,1);
     for i=1:nlinfams
        Ilin(i)=find(strcmp(All_fams,lin_fams{i}));
     end
     Ichroms = zeros(2,1);
     for i=1:2
        Ichroms(i)=find(strcmp(All_fams,chrom_fams{i}));
     end
   end
 end

if (not(strcmpi(optMode,'UC')))
  if (isfield(LatticeOptData,'All_fams'))

    Kall = zeros(LatticeOptData.nallfams,1);
    for i=1:LatticeOptData.nvars
        Kall(Idvs(i))=Ks(i);
    end
    for i=1:nlinfams
        Kall(Ilin(i))=MOGAResults.DVLins(i);
    end

    if (not(fitchromf)&&not(strcmpi(optMode,'UC')))
        for i=1:2
            Kall(Ichroms(i))=MOGAResults.ParetoFront(i0,nvars+nfitvars+6+i);
        end
    end

    Kall = getAllfams(2,ACHRO)
  end
end

  if (isfield(LatticeOptData,'All_famsO'))
    KallO = zeros(LatticeOptData.nallfamsO,1);
    for i=1:LatticeOptData.nvars
        KallO(Idvs(i))=Ks(i);
    end
    for i=1:nlinfams
        KallO(Ilin(i))=MOGAResults.DVLins(i);
    end
    if (not(fitchromf)&&not(strcmpi(optMode,'UC')))
        for i=1:2
            KallO(Ichroms(i))=MOGAResults.ParetoFront(i0,nvars+nfitvars+6+i);
        end
    end
  end
end

Knew = Ks;
if (not(strcmpi(optMode,'UC')))
    Kall_new = Kall;
    KallO_new = KallO;
else
    Kall_new=[];
    KallO_new=[];
    Kall = [];
end
%}
if (verboselevel>0)
    fprintf('**** \n');
    fprintf('%s Extracting individual n.%3d from MOGA set: %s \n', ...
             datetime, index, RunNumber);
    fprintf('**** \n');
end
%
% Sets the decision variables for the retrieved individual
%

UC   = setDVs(3, UC, LatticeOptData, Ks);
if (not(strcmpi(optMode,'UC')))
    LAT  = setDVs(2, LAT,LatticeOptData, Ks);
    if (isfield(LatticeOptData,'ACHROGRD'))
        ACHROGRD=setDVs(7,ACHROGRD,LatticeOptData, Ks);
    end
    Kall  = getAllfams(2,LAT,LatticeOptData);
    KallO = getAllfamsO(2,LAT,LatticeOptData);
else
    LAT  = UC;
    Kall = getAllfams(3,LAT,LatticeOptData);
    KallO = getAllfamsO(3,LAT,LatticeOptData);
end

%
% Loads linear lattice parameters for optimization mode 'SEXT'
%
if(strcmp(optMode,'SEXT')||strcmp(optMode,'NonLinear'))
  DVLins=MOGAResults.DVLins;
  LAT = setLins(2,LAT,LatticeOptData,DVLins);
  UC  = setLins(3,UC,LatticeOptData,DVLins);
  if (isfield(LatticeOptData,'ACHROGRD'))
    ACHROGRD=setLins(7,ACHROGRD,LatticeOptData, DVLins);
  end
  Kall  = getAllfams(2,LAT,LatticeOptData);
  KallO = getAllfamsO(2,LAT,LatticeOptData);
end

Kall_new = Kall;
KallO_new = KallO;

UC_tune=UC;
RunTitle = strcat(strrep(MOGAResults.RunNumber,'_','-'), ' Ind: ', num2str(index));


%% Lattice parameters/functions as derived from MOGA results directly
%
if (plotf)
    try
        figure;atplot(LAT);title(strcat('bare: ',RunTitle));
    catch ME
        fprintf('Error in atplot for MOGA lattice \n');
        fprintf('Error message was:%s \n',ME.message);
    end
end

try
    rpar=atsummary(LAT);
catch ME
    fprintf('%s Error in ExMOGA: atsummary full period from MOGA \n', datetime);
    fprintf('Error message was:%s \n',ME.message);
end

%% Unit cell
%
try
    rparaUC  = atsummary(UC);
    tunesuc0 = rparaUC.tunes;
    tunesuc1 = tunesuc0;
    [~, ~, ~, a2d, ~, ~] = Survey2D(UC,1.5*pi/180);
    totUCAngle = (a2d(1)-a2d(end))*180/pi;

catch ME
    fprintf('%s Error in ExMOGA: unit cell tune calculation for MOGA lattice\n', datetime);
    fprintf('Error message was:%s \n',ME.message);
    tunesuc0 = [NaN NaN];
end

%% Design Orbit
%
if (verboselevel>0)
    fprintf('%s ExMOGA: calculating design orbit \n', datetime);
end
if (split>1)
    LAT_SP     = splitlat(LAT,split);
else
    LAT_SP     = LAT;
end
[~, ~, ~, a2d, ~, ~] = Survey2D(LAT,0);
totAngle = (a2d(1)-a2d(end))*180/pi;

[s2d, x2d, y2d, a2d, baa, ban] = Survey2D(LAT_SP,totAngle/2*pi/180);
    
npo=numel(x2d);

geometry.DesignOrbit.s2d=s2d;
geometry.DesignOrbit.x2d=x2d;
geometry.DesignOrbit.y2d=y2d;
geometry.DesignOrbit.a2d=a2d;
geometry.DesignOrbit.baa=baa;
geometry.DesignOrbit.ban=ban;

if (not(isempty(ACHRO_ref)))
    if (split>1)
        ACHRO_ref_SP=splitlat(ACHRO_ref,split);
    else
        ACHRO_ref_SP=ACHRO_ref;
    end
    
    [s2d_ref, x2d_ref, y2d_ref, a2d_ref,~,~] = Survey2D(ACHRO_ref_SP,totAngle/2*pi/180);
    npo_ref = numel(x2d_ref);

    geometry.ref.DesignOrbit.s2d = s2d_ref;
    geometry.ref.DesignOrbit.x2d = x2d_ref;
    geometry.ref.DesignOrbit.y2d = y2d_ref;
    geometry.ref.DesignOrbit.a2d = a2d_ref;
    [~, ia, ~] = unique(x2d_ref);
    x2d_refu = x2d_ref(ia);
    y2d_refu = y2d_ref(ia);
    dist  = zeros(npo,1);
    for i = 1:npo
        yinter=interp1(x2d_refu,y2d_refu,x2d(i));
        dist(i)    = abs(yinter - y2d(i));
    end
else
    dist=nan(length(x2d),1);
    geometry.ref.DesignOrbit.s2d=nan(length(x2d),1);
    geometry.ref.DesignOrbit.x2d=nan(length(x2d),1);
    geometry.ref.DesignOrbit.y2d=nan(length(x2d),1);
    geometry.ref.DesignOrbit.a2d=nan(length(x2d),1); 
end

geometry.DesignOrbit.Deviation = dist;
LS.Lattice_Name=RunTitle;
LS.LattData.geometry = geometry;
LS.ACHROMAT = LAT;

if(plotf)
    plotGO(LS,'do','dodev');
end
%
% Field and gradient profiles
%
FG=calcFields(LAT,LatticeOptData.All_famsO,'split',split,'desc',RunTitle);
if (plotf)
    plotfield(FG);
end
LAT_tune=LAT;

%% Matching/Post processing
if (not(strcmpi(optMode,'UC')))
%% Matches unit cell tunes
LAT_UC=LAT;
Kold = Ks;
if (isfield(LatticeOptData,'All_fams'))
    Kall_old = Kall;
    Kall_new = Kall_old;
end
if (isfield(LatticeOptData,'All_famsO'))
    KallO_old = KallO;
    KallO_new = KallO_old;
end
%
%
if(fitucf&&not(isnan(tunesuc0(1)))&&not(isnan(tunesuc0(2))))
    if (verboselevel>0)
        fprintf('Fitting unit cell tunes from %4.3f %4.3f to %4.3f %4.3f \n',...
             tunesuc0(1), tunesuc0(2), tunesuc(1), tunesuc(2) );
    end
[UC_tune, its, penaltyuctune, ftunes]= fittuneRS(UC, tunesuc, uctune_fams{1}, uctune_fams{2}, 'maxits', nittune, 'Tol', TolTune, 'UseIntegerPart', false);
    try
        rpuc = atsummary(UC_tune);
        tunesuc1  = rpuc.tunes;
        if (verboselevel>0)
            fprintf('Unit cell tunes fit to[ %8.5f %8.5f ]\n',...
                 tunesuc1(1), tunesuc1(2));
             fprintf (' in %5d iterations and penalty = %8.2e \n', ...
                 its, penaltyuctune );
        end
        Knew     = getDVs(3, UC_tune, LatticeOptData); 
        LAT_UC   = setDVs(2, LAT, LatticeOptData, Knew);
 %
 % Here Knew contains NaNs as  the unit cell does not have all families
 % that are decision variables. 
 %
        Knew    = getDVs(2, LAT_UC, LatticeOptData);
        if (isfield(LatticeOptData,'All_fams'))
            Kall_new = getAllfams(2,LAT_UC,LatticeOptData);
            LAT_UC   = setAllfams(2,LAT_UC,LatticeOptData,Kall_new);
        end
        if isfield(LatticeOptData.ACHROGRD)
            ACHROGRD  = setDVs(7,ACHROGRD, LatticeOptData, Knew);
            KallO_new = getAllfamsO(7,ACHROGRD,LatticeOptData);
        end
 
    catch ME
        fprintf('Error in ExMOGA:atsummary of unit cell after unit cell fit \n');
        fprintf('Error message was:%s \n',ME.message);
        Dialog = questdlg('Revert lattice ?','Yes', 'No');
        switch Dialog
            case 'Yes'
                UC_tune=UC;
                LAT_UC = LAT;
                tunesuc1=tunesuc0;
                Knew = Kold;
                if (isfield(LatticeOptData,'All_fams'))
                    Kall_new=Kall_old;
                end
                if (isfield(LatticeOptData,'All_famsO'))
                    KallO_new = KallO_old;
                    ACHROGRD  = setAllfamsO(7,ACHROGRD,LatticeOptData,Kall_new);
                end
            case 'No'  
        end
    end
    %
    % Checks lattice of the full period with the  updated unit cell
    %
    try
        rpar=atsummary_fast(LAT_UC,isdipole);
    catch ME
        fprintf('Error in atsummary of full period after unit cell fit \n');
        fprintf('Error message was:%s \n',ME.message);
        Dialog = questdlg('Revert lattice ?','Yes', 'No');
        switch Dialog
            case 'Yes'
                UC_tune=UC;
                LAT_UC = LAT;
                tunesuc1=tunesuc0;
                Knew = Kold;
                if (isfield(LatticeOptData,'All_fams'))
                    Kall_new=Kall_old;
                end
                if (isfield(LatticeOptData,'All_famsO'))
                    KallO_new=KallO_old;
                    ACHROGRD=setAllfamsO(7,ACHROGRD,LatticeOptData,Kall_new);
                end
                fprintf('Lattice reverted \n');
            case 'No'
        end 
    end
end

%
%% Matches dispersion
%
Kold=Knew;

if (isfield(LatticeOptData,'All_fams'))
    Kall_old=Kall_new;
end
if (isfield(LatticeOptData,'All_famsO'))
    KallO_old=KallO_new;
end
LAT_md = LAT_UC;
   
if (fitdispf)
    if (verboselevel>0)
        fprintf('Matching dispersion \n');
    end
    Variab1 = atVariableBuilder(LAT_UC,...
               disp_fams{1},...
               {'PolynomB',{1,2}}); %TBC - add limits based on limits given to MOGA
%                lb(6),ub(6));
%    Variab2 = atVariableBuilder(LAT_UC,...
%               'reversebendmc2_sim',...
%               {'PolynomB',{1,2}},...
%               lb(7),ub(7));
%           
%    Variables = [Variab1;Variab2];
    Variables = Variab1;
    Constr1  = atlinconstraint(1,{{'Dispersion',{1}}},0,0,1);
  
    [LAT_md, penalty_md, dmin_md] = atmatch(LAT_UC,Variables,Constr1,...
                         TolDisp,nitdisp,0,@fminsearch); %  @fmincon  
    
    if (verboselevel>0)
        fprintf('Dispersion matched with penalty = %6.2e \n', sum(penalty_md.^2));
    end

    % Makes sure K and PolynomB are the same to avoid warning messages
    % later
    
    Knew   = getDVs(2,LAT_md,LatticeOptData);
    LAT_md = setDVs(2,LAT_md,LatticeOptData,Knew);

    if (isfield(LatticeOptData,'All_fams'))
        Kall_new = getAllfams(2,LAT_md,LatticeOptData);
        LAT_md   = setAllfams(2,LAT_md,LatticeOptData,Kall_new);
    end
    if (isfield(LatticeOptData,'All_famsO'))
        KallO_new = getAllfamsO(2,LAT_md,LatticeOptData);
        ACHROGRD  = setAllfamsO(7,ACHROGRD,LatticeOptData,KallO_new);
        KallO_new = getAllfamsO(7,ACHROGRD,LatticeOptData);
    end

    try
        rpar=atsummary_fast(LAT_md,isdipole);
    catch ME
        fprintf('Error in atsummary of full period after dispersion fit \n');
        fprintf('Error message was:%s \n',ME.message);
        Dialog = questdlg('Revert lattice ?','Yes', 'No');
        switch Dialog
            case 'Yes'
                LAT_md = LAT_UC;
                Knew = Kold;
                if (isfield(LatticeOptData,'All_fams'))
                    Kall_new=Kall_old;
                end
                if (isfield(LatticeOptData,'All_famsO'))
                    Kall_new=Kall_old;
                    ACHROGRD  = setAllfamsO(7,ACHROGRD,LatticeOptData,KallO_new);
                end
                fprintf('Lattice reverted \n');
                
            case 'No'
        end
    end
end

%% Matches beta functions at the centre of the long straight
%
Kold = Knew;
if (isfield(LatticeOptData,'All_fams'))
    Kall_old=Kall_new;
end
if (isfield(LatticeOptData,'All_famsO'))
    KallO_old=KallO_new;
end
LAT_bet = LAT_md;

if(fitbeta0f)
    if (verboselevel>0)
        fprintf('Fitting beta functions \n');
    end
    Variab1 = atVariableBuilder(LAT_bet,betaxy0_fams,...
                           {{'PolynomB',{1,2}},{'PolynomB',{1,2}}});% ,... TBC - add limits based on limits given to MOGA
%                           {lb(1),lb(2)},{ub(1),ub(2)});
    Constr1 = atlinconstraint(1,{{'beta',{1}},{'beta',{2}}},...
                  betas0,betas0,[1,1]);
    
    [LAT_bet, penalty_bet, dmin_bet] = atmatch(LAT_md,Variab1,Constr1,...
                         1E-6,1000,0,@fminsearch); %  @fmincon  
    if (verboselevel>0)   
        fprintf('Betas matched with penalty = %6.2e \n', sum(penalty_bet.^2));
    end
    
    Knew = getDVs(2,LAT_bet,LatticeOptData);
    LAT_bet = setDVs(2,LAT_bet,LatticeOptData,Knew);

    if (isfield(LatticeOptData,'All_fams'))
        Kall_new = getAllfams(2,LAT_bet,LatticeOptData);
        LAT_bet  = setAllfams(2,LAT_bet,LatticeOptData,Kall_new);
    end

    if (isfield(LatticeOptData,'All_famsO'))
        KallO_new = getAllfamsO(2,LAT_bet,LatticeOptData);
        ACHROGRD  = setAllfamsO(7,ACHROGRD,LatticeOptData,KallO_new);
        KallO_new = getAllfamsO(7,ACHROGRD,LatticeOptData);
    end
    
    try
        rpar=atsummary_fast(LAT_bet,isdipole);
    catch ME
        fprintf('Error in atsummary of full period after beta function fit \n');
        fprintf('Error message was: %s \n',ME.message);
        Dialog = questdlg('Revert lattice ?','Yes', 'No');
        switch Dialog
            case 'Yes'
                LAT_bet = LAT_md;
                Knew = Kold;
                 if (isfield(LatticeOptData,'All_fams'))
                    Kall_new=Kall_old;
                 end
                 if (isfield(LatticeOptData,'All_famsO'))
                    KallO_new=KallO_old;
                    ACHROGRD  = setAllfamsO(7,ACHROGRD,LatticeOptData,KallO_new);
                    KallO_new = getAllfamsO(7,ACHROGRD,LatticeOptData);
                end
                fprintf('Lattice reverted \n');
            case 'No'    
        end  
    end    
end
%
%% Matches tunes of the complete ring
%
Kold = Knew;
if (isfield(LatticeOptData,'All_fams'))
    Kall_old=Kall_new;
end
if (isfield(LatticeOptData,'All_famsO'))
    KallO_old=KallO_new;
end
LAT_tune = LAT_bet;

if(fittunef)
     qxfit = tunes(1)/Nper;
     qyfit = tunes(2)/Nper;
     if (verboselevel>0)
        fprintf('Fitting period tunes from [ %5.3f %5.3f ] to [ %5.3f %5.3f ] \n',...
              rpar.Qx_ring/Nper, rpar.Qy_ring/Nper,qxfit,qyfit );
     end
    
    [LAT_tune, its, penalty_tune]= fittuneRS(LAT_bet, [qxfit qyfit],...
         ringtune_fams{1}, ringtune_fams{2}, 'maxits', LatticeOptData.nittune, ...
         'Tol', LatticeOptData.TolTune,'UseIntegerPart',true);
    if (verboselevel>0)
        fprintf('Period Tune fit complete with penalty = %6.2e after %5d iterations \n', penalty_tune, its);
    end
 
    Knew = getDVs(2,LAT_tune,LatticeOptData);
    LAT_tune = setDVs(2,LAT_tune,LatticeOptData,Knew);

    if (isfield(LatticeOptData,'All_fams'))
        Kall_new = getAllfams(2,LAT_tune,LatticeOptData);
        LAT_tune = setAllfams(2,LAT_tune,LatticeOptData,Kall_new);
    end

    if (isfield(LatticeOptData,'All_famsO'))
        KallO_new = getAllfamsO(2,LAT_tune,LatticeOptData);
        ACHROGRD  = setAllfamsO(7,ACHROGRD,LatticeOptData,KallO_new);
        KallO_new = getAllfamsO(7,ACHROGRD,LatticeOptData);
    end
   
    try
        rpar=atsummary_fast(LAT_tune,isdipole);
        if (verboselevel>0)
            fprintf('Final ring tunes = [ %8.5f %8.5f ] \n', rpar.Qx_ring, rpar.Qy_ring);
        end
    
     catch ME
        fprintf('Error in atsummary of full period after tune fit \n');
        fprintf('Error message was: %s \n',ME.message);
        Dialog = questdlg('Revert lattice ?','Yes', 'No');
        switch Dialog
            case 'Yes'
                LAT_tune = LAT_bet;
                Knew = Kold;
                if (isfield(LatticeOptData,'All_fams'))
                    Kall_new=Kall_old;
                end
                if (isfield(LatticeOptData,'All_famsO'))
                    Kall_new=Kall_old;
                    ACHROGRD  = setAllfamsO(7,ACHROGRD,LatticeOptData,KallO_new);
                end
                fprintf('Lattice reverted \n');
                
            case 'No'
        end 
    end    
end

%% Plots matched Lattice Functions
if (plotf)
    try
      figure;atplot(LAT_tune);title(strcat('matched: ',RunTitle)); 
    catch ME
        fprintf('Error in atplot of final lattice \n');
        fprintf('Error message was: %s \n',ME.message);
    end
end
end
%
%% Matches chromaticity 
%
if (isfield(LatticeOptData,'All_fams'))
    Kall_old=Kall_new;
end
if (isfield(LatticeOptData,'All_famsO'))
    KallO_old=KallO_new;
end

if(fitchromf&&~isempty(chrom_fams))
     if (verboselevel>0) 
        fprintf('%s ExMOGA: matching chromaticity... \n', datetime); 
     end
     try 
        [LAT_C, Penalty, its]=fitchroit(LAT_tune, chrom_fams, chroms0, Nitchro, TolChrom); 
        I_sc1  = find(atgetcells(LAT_C, 'FamName', chrom_fams{1}));
        K_sc1  = atgetfieldvalues(LAT_C, I_sc1, 'PolynomB', {3});
        I_sc2  = find(atgetcells(LAT_C, 'FamName', chrom_fams{2}));
        K_sc2  = atgetfieldvalues(LAT_C, I_sc2, 'PolynomB', {3});
        Sc1    = K_sc1(1);
        Sc2    = K_sc2(1);
        if (verboselevel>0)
            fprintf('%s ExMOGA: chromaticity matched to %5.2 / %5.2f with penalty = %6.2e in %2d iterations\n', ...
                datetime, chroms0(1), chroms0(2), Penalty, its);
        end

        if(strcmpi(optMode,'UC'))
            nLAT=3;
        else
            nLAT=2;
        end

        if (isfield(LatticeOptData,'All_fams'))
            Kall_new = getAllfams(nLAT,LAT_C,LatticeOptData);
            LAT_C    = setAllfams(nLAT,LAT_C,LatticeOptData,Kall_new);
        end
        if (isfield(LatticeOptData,'All_famsO'))
            KallO_new = getAllfamsO(nLAT,LAT_C,LatticeOptData);
            if (not(isempty(ACHROGRD)))
                ACHROGRD  = setAllfamsO(7,ACHROGRD,LatticeOptData,KallO_new);
                KallO_new = getAllfamsO(7,ACHROGRD,LatticeOptData);
            end
        end
     catch ME
        fprintf('%s ExMOGA: error in chromaticity fit \n', datetime);
        fprintf('Error message was: %s \n',ME.message);
        LAT_C = LAT_tune;
        if (not(isempty(ACHORGRD)))
            ACHROGRD  = setAllfamsO(7,ACHROGRD,LatticeOptData,KallO_old);
        end
     end
     LAT_tune=LAT_C;
     
     try
        rpar=atsummary(LAT_tune);
     catch ME
        fprintf('Error in atsummary of full period after chromaticity fit \n');
        fprintf('Error message was: %s \n',ME.message);
        Sc1=NaN;
        Sc2=NaN;
    end    
else
    Sc1=NaN;
    Sc2=NaN;
end

%% Calculates and Plots the Dynamic Aperture
%
if (plotdaf)
    DAS=calcDA(ACHROGRD,DAoptions,'plot','4d');
    fprintf('Dynamic Aperture = %4.2f mm**2 \n',DAS.outputs.DA);
else
    DAS=[];
end

%
%% Full ring tune scan - done on the lattice without octupoles
% 
if (tunescanf)
    LAT_scan=LAT_tune;
    Qxmin = Qrange(1);
    Qxmax = Qrange(2);
    Qymin = Qrange(3);
    Qymax = Qrange(4);

    qxmin = qrange(1);
    qxmax = qrange(2);
    qymin = qrange(3);
    qymax = qrange(4);

    Npqx = Npq(1);
    Npqy = Npq(2);

    nx = (Qxmax-Qxmin+1)*Npqx;
    ny = (Qymax-Qymin+1)*Npqy;
    dqx = (qxmax-qxmin)/(Npqx-1);
    dqy = (qymax-qymin)/(Npqy-1);

    ntot = nx*ny;
    fprintf('Starting tune scan with %3d grid points \n', ntot);
    fb=waitbar(0,'Starting Tune Scan...', 'Name','Tune Scan Progress', 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(fb,'canceling',0);
    I_sc1 = find(atgetcells(LAT_tune, 'FamName', chrom_fams{1}));
    I_sc2 = find(atgetcells(LAT_tune, 'FamName', chrom_fams{2}));

 
    DAScan = zeros(nx,ny);
    EmitScan = zeros(nx,ny);
    BetaX0Scan = zeros(nx,ny);
    BetaY0Scan = zeros(nx,ny);
    JxScan     = zeros(nx,ny);
    qxscan     = zeros(1,nx);
    qyscan     = zeros(1,ny);

    k=1;
    for i=Qxmin:Qxmax
        for j=1:Npqx
            qxscan(k)=i+qxmin+dqx*(j-1);
            k=k+1;
        end
    end

    k=1;
    for i=Qymin:Qymax
        for j=1:Npqy
            qyscan(k)=i+qymin+dqy*(j-1);
            k=k+1;
        end
    end

    
    nfrac=0;
    nscan=0;

    for i=1:nx
        Qxfit = qxscan(i);
        qxfit = Qxfit/Nper;
        for j=1:ny
            Qyfit = qyscan(j);
            qyfit=Qyfit/Nper;
            try
                [LAT_scan, its, penalty_tune]= fittuneRS(LAT_scan, [qxfit qyfit],...
                                        ringtune_fams{1}, ringtune_fams{2},...
                                        'maxits', LatticeOptData.nittune,...
                                        'Tol', LatticeOptData.TolTune,'UseIntegerPart', true);
                if (verbose)
                    fprintf('Ring Tunes fit to [ %4.2f , %4.2f ] with penalty = %6.2e after %5d iterations \n', Qxfit, Qyfit, penalty_tune, its);
                end
            catch ME
                fprintf('ExMOGA Error in tune fit at ring tunes  [ %4.2f , %4.2f ] \n', Qxfit, Qyfit);
                fprintf('Error message was:%s \n',ME.message);
                LAT_scan=LAT_tune;
                EmitScan(i,j)   = NaN;
                BetaX0Scan(i,j) = NaN;
                BetaY0Scan(i,j) = NaN;
                JxScan(i,j)     = NaN;
                DAScan(i,j)     = NaN;
                nscan=nscan+1;
                nfracnew=nscan/ntot*100;
                if ((nfracnew-nfrac)>1)
                    waitbar(nscan/ntot,fb,strcat(sprintf('%3.0f %s',nfracnew,'%')));
                    nfrac=nfracnew;
                end
                continue
            end
            try
                rparascan=atsummary_fast(LAT_scan,isdipole);
                EmitScan(i,j)   = rparascan.naturalEmittance;
                BetaX0Scan(i,j) = rparascan.beta0(1);
                BetaY0Scan(i,j) = rparascan.beta0(2);
                JxScan(i,j)     = rparascan.damping(1);
            catch ME
                fprintf('ExMOGA Error in atsummary at ring tunes  [ %4.2f , %4.2f ] \n', Qxfit, Qyfit);
                fprintf('Error message was:%s \n',ME.message);
                EmitScan(i,j)   = NaN;
                BetaX0Scan(i,j) = NaN;
                BetaY0Scan(i,j) = NaN;
                JxScan(i,j)     = NaN;
                DAScan(i,j)     = NaN;
                nscan=nscan+1;
                nfracnew=nscan/ntot*100;
                if ((nfracnew-nfrac)>1)
                    waitbar(nscan/ntot,fb,strcat(sprintf('%3.0f %s',nfracnew,'%')));
                    nfrac=nfracnew;
                end
                continue
            end
            %
            % Fits chromaticity 
            %
            try 
                [LAT_scan, Penalty, its]=fitchroit(LAT_scan, chrom_fams, chroms0, Nitchro, TolChrom); 
                K_sc1  = atgetfieldvalues(LAT_scan, I_sc1, 'PolynomB', {3});
                K_sc2  = atgetfieldvalues(LAT_scan, I_sc2, 'PolynomB', {3});
                Sc1    = K_sc1(1);
                Sc2    = K_sc2(1);
            if (verboselevel>0)
                fprintf('Chromaticity matched with penalty = %6.2e in %2d iterations\n', Penalty, its);
            end
            catch ME
                fprintf('Error in ExMOGA: chromaticity fit \n');
                fprintf('Error message was: %s \n',ME.message);
                LAT_scan = LAT_tune;
                Sc1=NaN;
                Sc2=NaN;
                DAScan(i,j) = NaN;
                nscan=nscan+1;
                nfracnew=nscan/ntot*100;
                if ((nfracnew-nfrac)>1)
                    waitbar(nscan/ntot,fb,strcat(sprintf('%3.0f %s',nfracnew,'%')));
                    nfrac=nfracnew;
                end
                continue
            end
              
            if (not(isnan(Sc1)&&not(isnan(Sc2))))
                [DA,~]=calcDA(LAT_scan,DAoptions);
            else
                DA=NaN;
            end
            DAScan(i,j)=DA;
            if (verboselevel>0)
                fprintf('i = %3d j = %3d qx= %5.2f qy= %5.2f DA = %5.2f mm**2 Emit = %5.2f pmrad \n', i, j, Qxfit, Qyfit, DA, EmitScan(i,j)*1e12);
            end
            nscan=nscan+1;
            nfracnew=nscan/ntot*100;
            if ((nfracnew-nfrac)>1)
                    waitbar(nscan/ntot,fb,strcat(sprintf('%3.0f %s',nfracnew,'%')));
                    nfrac=nfracnew;
            end
            if getappdata(fb,'canceling')
                break
            end
        end
        if getappdata(fb,'canceling')
           fprintf ('Canceling at nscan = %5d i = %3d j = %3d qx= %5.2f qy= %5.2f \n', nscan, i, j, Qxfit, Qyfit);
           break
        end
    end 
    delete(fb);
    rp.outputs.qxscan=qxscan;
    rp.outputs.qyscan=qyscan;
    rp.outputs.DAScan=DAScan;
    rp.outputs.EmitScan=EmitScan*1e12;
    rp.outputs.BetaX0Scan = BetaX0Scan;
    rp.outputs.BetaY0Scan = BetaY0Scan;
    rp.outputs.JxScan = JxScan;
    if (strcmp(plottunescanf,'Y'))
        PlotEMTuneScan(rp);
    end
    [mdaj,jmaxda] = max(DAScan');[maxda,imaxda] = max(mdaj);
    rp.outputs.DAmax = maxda;
    rp.outputs.qxmaxDA = qxscan(imaxda);
    rp.outputs.qymaxDA = qyscan(jmaxda(imaxda));
end


%% Collect data in output structure
%
% "inputs" field documents input parameters to ExMOGA
% "outputs" field documents output parameters 
%

rp.inputs.MOGAResults=MOGAResults;
rp.inputs.MOGAIndex=index;
rp.inputs.fitucf= fitucf; 
rp.inputs.fitdispf  = fitdispf; 
rp.inputs.fitbeta0f = fitbeta0f; 
rp.inputs.fittunef  = fittunef;
rp.inputs.fitchromf = fitchromf;
rp.inputs.tunescanf = tunescanf;
rp.inputs.saveOPAf  = saveOPAf;
rp.inputs.tunesuc = tunesuc;
rp.inputs.betas0=betas0;
rp.inputs.tunes=tunes;
rp.inputs.Qrange=Qrange;
rp.inputs.qrange=qrange;
rp.inputs.Npq=Npq;
rp.inputs.DAoptions=DAoptions;

rp.outputs.atsummary = rpar;
rp.outputs.K0        = Ks;
rp.outputs.Knew      = Knew;
rp.outputs.tunesuc0  = tunesuc0;
rp.outputs.tunesuc1  = tunesuc1;
rp.outputs.index     = index;
rp.outputs.DAoptions = DAoptions;
rp.outputs.LatticeOptData = LatticeOptData;
rp.outputs.Sc1 = Sc1;
rp.outputs.Sc2 = Sc2;
rp.outputs.totUCAngle   = totUCAngle;
rp.outputs.totAngle = totAngle;
if (not(isempty(DAS)))
    rp.outputs.DA=DAS.outputs.DA;
else
    rp.outputs.DA=nan;
end

if(strcmp(optMode,'SEXT')||strcmp(optMode,'NonLinear'))
    rp.outputs.DVLins = DVLins;
end

if(isfield(LatticeOptData,'All_fams'))
    rp.outputs.Kall = Kall_new;
end

if(isfield(LatticeOptData,'All_famsO'))
    rp.outputs.KallO = KallO_new;
end

if (not(strcmpi(optMode,'UC')))
    rp.outputs.ACHRO = LAT_tune;
end

rp.outputs.UC=UC_tune;

if (not(strcmpi(optMode,'UC')))
 if (isfield(LatticeOptData,'RINGGRD'))
    if (isfield(LatticeOptData,'All_famsO'))
        rp.outputs.RINGGRD=setAllfamsO(6,RINGGRD,LatticeOptData,KallO_new);
    else
        rp.outputs.RINGGRD=setAllfams(6,RINGGRD,LatticeOptData,Kall_new);
    end
 end

 if (isfield(LatticeOptData,'ACHROGRD'))
    if (isfield(LatticeOptData,'All_famsO'))
        rp.outputs.ACHROGRD = setAllfamsO(7,ACHROGRD,LatticeOptData,KallO_new);
    else
        rp.outputs.ACHROGRD = setAllfams(7,ACHROGRD,LatticeOptData,Kall_new);
    end
 end
end



%
%% Saves OPA file
%
if(saveOPAf)
    if (isfield(LatticeOptData,'All_famsO'))
        saveOPA(KallO_new,[],[],LatticeOptData,RunNumber,'oct');
    else
        if (isfield(LatticeOptData,'All_fams'))
            saveOPA(Kall_new,[],[],LatticeOptData,RunNumber);
        else
            saveOPA([],Knew,[Sc1 Sc2],LatticeOptData,RunNumber);
        end
    end
end

if ((verboselevel>0)||tunescanf)
    toc;
end
%% Auxiliary functions
function fp=FractionalPart(x)
    fp=x-fix(x);




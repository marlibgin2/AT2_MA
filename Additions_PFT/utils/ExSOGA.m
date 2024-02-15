function rpara=ExSOGA(varargin)
% 
% Examines result of SOGA Scan, produces optics function plots 
%
%% inputs:
% mandatory arguments:
%           SOGAResults: structure containig SOGA results (see SOGA.m
%           index     : index of line in "Population" an "Scores" fields in 
%                       SOGAResults structure from which data will be 
%                       retrieved. If negative, the line of
%                       lowest objective function is chosen. If zero, the
%                       values in field OptDVs of the SOGAResults
%                       structure are chosen.
%
% Optional arguments: presence in the argument list activates the option
%           plot     : if "Y", plots lattice functions for the extracted lattice
%                       as well for the descendent version with dispersion and tune fits
%           fitucf    : if 'Y' fits the unit cell tunes to specified values
%           fitdispf  : if 'Y' matches dispersion to zero in long straight
%           fitbeta0f : if 'Y' fits the horizontal and vertical beta
%                       functions at centre of long stright
%           fittunef  : if "Y" fits lattice to match dispersion and fit period tunes
%                       all other fits.
%           fitchromf : of "Y fist chromatrocyt with families listed in the chrom_fams
%           plotdaf   : if 'Y plots calculates and plots the on-energy
%                       dynamic  aperture.
%           betas0    : 1x2 matrix of desired long straight hor/vert beta
%                       functions.
%           tunesuc   : 1x2 matrix of desired unit cell fractional tunes. (isolated cell)    
%           tunes     : 1x2 matrix of desired full ring fractional tunes.
%           saveOPAf  : if 'Y' saves a text file with the decision variables set for OPA run. 
%                       the file is based on a template - There atre two
%                       templates, one used in "simplified" mode (dipoles hae a single slice)and the other
%                       in "complete" mode (dipoles are sliced as in the
%                       standard lattice.
%
% outputs:
%          rpara: structure containing lattice properties as produced by
%                 function atsummary_fast.m
%
%
%
%% Input argument parsing
[SOGAResults,index]=getargs(varargin,[],1);
plotf          = any(strcmpi(varargin,'plot'));
fitucf         = any(strcmpi(varargin,'fituc'));
fitdispf       = any(strcmpi(varargin,'fitdisp')); 
fitbeta0f      = any(strcmpi(varargin,'fitbeta0')); 
fittunef       = any(strcmpi(varargin,'fittune'));
fitchromf      = any(strcmpi(varargin,'fitchrom')); 
plotdaf        = any(strcmpi(varargin,'plotda')); 
tunescanf      = any(strcmpi(varargin,'tunescan'));
plottunescanf  = any(strcmpi(varargin,'plottunscan'));
saveOPAf       = any(strcmpi(varargin,'saveOPA'));
verbose        = any(strcmpi(varargin,'verbose'));

tunesuc        = getoption(varargin,'tunesuc',[3/7 1/7]);
betas0         = getoption(varargin,'betas0',[9.0 2.0]);
tunes          = getoption(varargin,'tunes',[46.20 16.28]);
DAoptions      = getoption(varargin,'DAoptions',SOGAResults.LatticeOptData.DAoptions);
%% Preamble
LatticeOptData = SOGAResults.LatticeOptData;

isdipole    = LatticeOptData.isdipole;
nvars       = LatticeOptData.nvars;
LAT         = LatticeOptData.ACHRO;
UC          = LatticeOptData.UC;
Nper        = LatticeOptData.Nper;
scan_fams   = LatticeOptData.scan_fams;

RunNumber   = SOGAResults.RunNumber;
lb          = SOGAResults.lb;
ub          = SOGAResults.ub;

Sc1 = nan;
Sc2 = nan;

%
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
    nfitvars    = size(SOGAResults.fitvars,1);
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

if (isfield(LatticeOptData,'All_fams'))
    All_fams = LatticeOptData.All_fams;
    Iuctune = zeros(size(uctune_fams,1));
    for i=1:size(uctune_fams,1)
        Iuctune(i)=find(strcmp(All_fams,uctune_fams{i}));
    end
    Idisp = zeros(size(disp_fams,1));
    for i=1:size(size(disp_fams,1))
        Idisp(i)=find(strcmp(All_fams,disp_fams{i}));
    end
    Ibetaxy0 = zeros(size(betaxy0_fams,1));
    for i=1:size(betaxy0_fams,1)
        Ibetaxy0(i) = find(strcmp(All_fams,betaxy0_fams{i}));
    end
    Iringtune = zeros(size(ringtune_fams,1));
    for i=1:size(ringtune_fams,1)
        Iringtune(i)=find(strcmp(All_fams,ringtune_fams{i}));
    end
    Idvs = zeros(LatticeOptData.nvars);
    for i=1:nvars
        Idvs(i)=find(strcmp(All_fams,scan_fams{i}));
    end
    Ilin = zeros(nlinfams);
    for i=1:nlinfams
         Ilin(i)=find(strcmp(All_fams,lin_fams{i}));
    end
    Ichroms = zeros(2);
    for i=1:2
         Ichroms(i)=find(strcmp(All_fams,chrom_fams{i}));
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


%% Parameters for dynamic aperture calculation
% note that these may be different from the parameters at the time
% SOGA was run and recorded in the structure 
% SOGAResults.LatticeOptData.DAoptions
% 
%
% If the argumetn DAOptions given to the funcion is empty, the
% corresponding DAOptions from 
%

chroms0      = DAoptions.chroms0; % Target chromaticity for one superperiod
TolChrom     = DAoptions.TolChrom;% Chromaticity tolerances
Nitchro      = DAoptions.Nitchro; % Max n. iterations of chromaticty correction


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

if (strcmp(DAoptions.DAmode,'grid')&&nargin>nargs)
  
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


%
%% Searches for desired individual in the final population
%
if (index==0)
        Ks = SOGAResults.optDVs;
        im = NaN;
else
    if (index<0)
        [~,im] = min(SOGAResults.scores);
    else
        im = index;
    end
    Ks = SOGAResults.Population(im,1:nvars);
end

Knew = Ks;

if (isfield(LatticeOptData,'All_fams'))
    Kall = zeros(LatticeOptData.nallfams,1);
    for i=1:LatticeOptData.nvars
        Kall(Idvs(i))=Ks(i);
    end
    for i=1:nlinfams
        Kall(Ilin(i))=SOGAResults.DVLins(i);
    end
%    I_famsAllF  = LatticeOptData.IfamsAllF;
%    I_famsAllUC = LatticeOptData.IfamsAllUC;
end
Kall_new = Kall;


if (verbose)
    fprintf('**** \n');
    if (index==0)
        fprintf('%s Extracting optimum individual from SOGA set: %s \n', datetime, RunNumber);
        fprintf('with objective function = %5.2f \n', SOGAResults.score);
    else
        fprintf('%s Extracting individual %d3 from SOGA set: %s \n', datetime, im, RunNumber);
        fprintf('with objective function = %5.2f \n', SOGAResults.scores(im));
    end

    fprintf('**** \n');
    tic;
end

%
%% Sets the decision variables for the retrieved individual
%

LAT  = setDVs(2, LAT,LatticeOptData, Ks);
UC   = setDVs(3, UC, LatticeOptData, Ks);

% loads linear lattice parameters for optimization mode 'SEXT' (legacy) or
% 'NonLinear'

if(strcmp(optMode,'NonLinear')||strcmp(optMode,'SEXT'))
  DVLins=SOGAResults.DVLins;
  LAT = setLins(2,LAT,LatticeOptData,DVLins);
% UC  = setLins(3,UC,LatticeOptData,DVLins);
end
%
%% Fits chromaticity 
%
Sc1    = nan;
Sc2    = nan;
LAT_C=LAT;
if(fitchromf)
     if (verbose) 
        fprintf('Fitting Chromaticity... \n'); 
     end
     
     try 
        [LAT_C, Penalty, its]=fitchroit(LAT, chrom_fams, chroms0, Nitchro, TolChrom); 
        I_sc1  = find(atgetcells(LAT_C, 'FamName', chrom_fams{1}));
        K_sc1  = atgetfieldvalues(LAT_C, I_sc1, 'PolynomB', {3});
        I_sc2  = find(atgetcells(LAT_C, 'FamName', chrom_fams{2}));
        K_sc2  = atgetfieldvalues(LAT_C, I_sc2, 'PolynomB', {3});
        Sc1    = K_sc1(1);
        Sc2    = K_sc2(1);
        if (isfield(LatticeOptData,'All_fams'))
            Kall_new(Ichroms(1))=Sc1;
            Kall_new(Ichroms(2))=Sc2;
        end
        if (verbose)
            fprintf('Chromaticty matched with penalty = %6.2e in %2d iterations\n', Penalty, its);
        end
     catch ME
        fprintf('Error in ExSOGA: chromaticity fit \n');
        fprintf('Error message was: %s \n',ME.message);
        LAT_C = LAT;
     end
end
%
%% Plots lattice functions as derived from SOGA results directly
%
if (plotf)
    
    try
        figure;atplot(LAT_C);title('SOGA bare');
  
    catch ME
          fprintf('Error atplot for SOGA lattice \n');
          fprintf('Error message was:%s \n',ME.message);
    end
end

try
    rpara=atsummary_fast(LAT_C,isdipole);
catch ME
    fprintf('Error in atsummary full period from SOGA \n');
    fprintf('Error message was:%s \n',ME.message);
end
 
%% legacy from ExMOGA
%LAT_tune = LAT;
%try
%    rparaUC = atsummary(UC);
%    tunesuc0 = rparaUC.tunes;
%catch
%    fprintf('error in unit cell tune calculation \n');
%    tunesuc0 = [NaN NaN];
%end

%tunesuc1=tunesuc0;

%LAT_UC=LAT;
%Kold = Ks;

%
%linear lattice fits are only performed if mode is not "SEXT"
%
%if(not(strcmp(optMode,'SEXT')))
% if(strcmp(fitucf,'Y')&&not(isnan(tunesuc0(1)))&&not(isnan(tunesuc0(2))))
%    if (strcmp(verbose,'Y'))
%        fprintf('Fitting unit cell tunes from %4.3f %4.3f to %4.3f %4.3f \n',...
%             tunesuc0(1), tunesuc0(2), tunesuc(1), tunesuc(2) );
%    end
%    fam1 = 'reversebend_sim';
%    if (strcmpi(optMode,'SIMP'))
%        fam2 = 'sshdip';
%    else
%        fam2 = 'dip';
%    end
%   UC_tune = atfittune(UC,  tunesuc, fam1, fam2);
%   UC_tune = atfittune(UC_tune,tunesuc , fam1, fam2);
%    [UC_tune, its, penaltyuctune ftunes]= fittuneRS(UC, tunesuc, fam1, fam2, 1000, 1E-7, 'N');
%   [UC_tune, penaltyuctune] = fitfullTune(UC, fam1, fam2, lb(3:4), ub(3:4), tunesuc);
%    try
%        rpuc = atsummary(UC_tune);
%        tunesuc1  = rpuc.tunes;
%        if (strcmpi(verbose,'Y'))
%            fprintf('Unit cell tunes fit to[ %8.5f %8.5f ]\n',...
%                 tunesuc1(1), tunesuc1(2));
%             fprintf (' in %5d iterations and penalty = %8.2e \n', ...
%                 its, penaltyuctune );
%        end
%        Knew    = getDVs(3, UC_tune, LatticeOptData); 
%        UC_tune = setDVs(3, UC_tune, LatticeOptData,Knew); 
%        LAT_UC  = setDVs(2, LAT, LatticeOptData, Knew);
%
 % Here Knew contains NaNs as  the unit cell does not have all families.
 % Below we fix this.
%        Knew    = getDVs(2, LAT_UC, LatticeOptData);
%                    
%    catch ME
%        fprintf('Error in atsummary of unit cell after unit cell fit \n');
%        fprintf('Error message was:%s \n',ME.message);
%        Dialog = questdlg('Revert lattice ?','Yes', 'No');
%        switch Dialog
%            case 'Yes'
%                UC_tune=UC;
%                LAT_UC = LAT;
%                tunesuc1=tunesuc0;
%                Knew = Kold;
%                
%            case 'No'
%                
%        end % switch
%    end
 
    %
    % Checks lattice of the full period with the  updated unit cell
    %
%    try
%        rpara=atsummary_fast(LAT_UC,isdipole);
%    catch ME
%        fprintf('Error in atsummary of full period after unit cell fit \n');
%        fprintf('Error message was:%s \n',ME.message);
%        Dialog = questdlg('Revert lattice ?','Yes', 'No');
%        switch Dialog
%            case 'Yes'
%                UC_tune=UC;
%                LAT_UC = LAT;
%                tunesuc1=tunesuc0;
%                Knew = Kold;
%                fprintf('Lattice reverted \n');
%            case 'No'
% 
%        end % switch    
%    end
%end

%Kold=Knew;
%LAT_md = LAT_UC;
   
%if (strcmp(fitdispf,'Y'))
%    if (strcmp(verbose,'Y'))
%        fprintf('Matching dispersion \n');
%    end
%    switch optMode
%        case {'SIMP','COMP'}
%            rbfam = 'reversebendmc1_sim';
%        case 'CHRO'
%            rbfam = 'reversebendmc1_chro';
%    end
%    Variab1 = atVariableBuilder(LAT_UC,...
%               rbfam,...
%               {'PolynomB',{1,2}},...
%                lb(6),ub(6));
% 
%    Variab2 = atVariableBuilder(LAT_UC,...
%               'reversebendmc2_sim',...
%               {'PolynomB',{1,2}},...
%               lb(7),ub(7));
%           
%    Variables = [Variab1;Variab2];
%    Variables = Variab1;
%    Constr1  = atlinconstraint(1,{{'Dispersion',{1}}},0,0,1);
  
%    [LAT_md, penalty_md, dmin_md] = atmatch(LAT_UC,Variables,Constr1,...
%                         1E-8,1000,0,@fminsearch); %  @fmincon  
    
%    if (strcmp(verbose,'Y'))
%        fprintf('Dispersion matched with penalty = %6.2e \n', sum(penalty_md.^2));
%    end

    % Makes sure K and PolynomB are the same to avoid warning messages
    % later
    
%    Knew = getDVs(2,LAT_md,LatticeOptData);
%    LAT_md = setDVs(2,LAT_md,LatticeOptData,Knew);
    
    
%    try
%        rpara=atsummary_fast(LAT_md,isdipole);
%    catch ME
%        fprintf('Error in atsummary of full period after dispersion fit \n');
%        fprintf('Error message was:%s \n',ME.message);
%        Dialog = questdlg('Revert lattice ?','Yes', 'No');
%        switch Dialog
%            case 'Yes'
%                LAT_md = LAT_UC;
%                Knew = Kold;
%                fprintf('Lattice reverted \n');
%                
%            case 'No'

%        end % switch 
%    end
%end

%Kold = Knew;
%LAT_bet = LAT_md;

%if(strcmp(fitbeta0f,'Y'))
%    if (strcmp(verbose,'Y'))
%        fprintf('Fitting beta functions \n');
%    end
%    Variab1 = atVariableBuilder(LAT_bet,{'qfend_sim','qdend_sim'},...
%                           {{'PolynomB',{1,2}},{'PolynomB',{1,2}}},...
%                           {lb(1),lb(2)},{ub(1),ub(2)});
%    Constr1 = atlinconstraint(1,{{'beta',{1}},{'beta',{2}}},...
%                  betas0,betas0,[1,1]);
%    
%    [LAT_bet, penalty_bet, dmin_bet] = atmatch(LAT_md,Variab1,Constr1,...
%                         1E-6,1000,0,@fminsearch); %  @fmincon  
%    if (strcmp(verbose,'Y'))   
%        fprintf('Betas matched with penalty = %6.2e \n', sum(penalty_bet.^2));
%    end
    
%    Knew = getDVs(2,LAT_bet,LatticeOptData);
%    LAT_bet = setDVs(2,LAT_bet,LatticeOptData,Knew);
    
%    try
%        rpara=atsummary_fast(LAT_bet,isdipole);
%    catch ME
%        fprintf('Error in atsummary of full period after beta function fit \n');
%        fprintf('Error message was: %s \n',ME.message);
%        Dialog = questdlg('Revert lattice ?','Yes', 'No');
%        switch Dialog
%            case 'Yes'
%                LAT_bet = LAT_md;
%                Knew = Kold;
%                fprintf('Lattice reverted \n');
%                
%            case 'No'
%                
%        end % switch 
%    end    
%end 
    
%Kold = Knew;
%LAT_tune = LAT_bet;
%if(strcmp(fittunef,'Y'))
%     qxfit = tunes(1)/Nper;
%     qyfit = tunes(2)/Nper;
%     if (strcmp(verbose,'Y'))
%        fprintf('Fitting period tunes from [ %5.3f %5.3f ] to [ %5.3f %5.3f ] \n',...
%              rpara.Qx_ring/Nper, rpara.Qy_ring/Nper,qxfit,qyfit );
%     end
%    
%     LAT_tune = atfittune(LAT_md,  [qxfit qyfit], 'qfend_sim', 'qdend_sim', 'UseIntegerPart');
%     LAT_tune = atfittune(LAT_tune,[qxfit qyfit] ,'qfend_sim', 'qdend_sim', 'UseIntegerPart');

     %lb = MOGAResults.lb(1:2);
     %ub = MOGAResults.ub(1:2);       
 %   [LAT_tune, penalty_tune] = fitfullTune(LAT_bet, 'qfend_sim', 'qdend_sim', lb, ub, [qxfit qyfit]);
 %    [LAT_tune, its, penalty_tune]= fittuneRS(LAT_bet, [qxfit qyfit],...
 %        'qfend_sim', 'qdend_sim', 1000, 1E-7,'Y');
 %   if (strcmpi(verbose,'Y'))
 %       fprintf('Period Tune fit complete with penalty = %6.2e after %5d iterations \n', penalty_tune, its);
 %   end
% 
%    Knew = getDVs(2,LAT_tune,LatticeOptData);
%    LAT_tune = setDVs(2,LAT_tune,LatticeOptData,Knew);
   
%    try
%        rpara=atsummary_fast(LAT_tune,isdipole);
%        if (strcmp(verbose,'Y'))
%            fprintf('Final ring tunes = [ %8.5f %8.5f ] \n', rpara.Qx_ring, rpara.Qy_ring);
%        end
    
%     catch ME
%       fprintf('Error in atsummary of full period after tune fit \n');
%        fprintf('Error message was: %s \n',ME.message);
%        Dialog = questdlg('Revert lattice ?','Yes', 'No');
%        switch Dialog
%            case 'Yes'
%                LAT_tune = LAT_bet;
%                Knew = Kold;
%                fprintf('Lattice reverted \n');
%                
%            case 'No'
                
%        end % switch 
%    end    
% end
%else
%    LAT_tune=LAT;
%end

%if (strcmp(plotf,'Y'))
%    
%    try
%        figure;atplot(LAT);title('MOGA bare');
%  
%    catch ME
%          fprintf('Error atlinopt for SOGA lattice \n');
%          fprintf('Error message was:%s \n',ME.message);
%    end
%end

%
%% Calculates and plots Dynamic Aperture
%
if (plotdaf)
    DA=CalcPlotDA(LAT_C,DAoptions,plotdaf);
    fprintf('Dynamic Aperture = %4.2f mm**2 \n',DA);
else
    DA=NaN;
end

%% Adds info to output structure

rpara.Knew  = Knew;
rpara.K0    = Ks;
rpara.Sc1   = Sc1;
rpara.Sc2   = Sc2;
rpara.Kall  = Kall_new;
rpara.DA    = DA;
%rpara.tunesuc0=tunesuc0;
%rpara.tunesuc1=tunesuc1;
%rpara.index=iindex;
rpara.LatticeOptData = LatticeOptData;
rpara.SOGARunNumber = RunNumber;
rpara.SOGAIndex=index;
rpara.SOGAIndividual=im;

if(strcmp(optMode,'NonLinear')||strcmp(optMode,'SEXT'))
    rpara.DVLins = DVLins;
end
rpara.ACHRO = LAT_C;

%
%% Saves OPA file
%
if(saveOPAf)
    if (isfield(LatticeOptData,'All_fams'))
        saveOPA(Kall_new,[],[],LatticeOptData,RunNumber);
    else
        saveOPA([],Knew,[Sc1 Sc2],LatticeOptData,RunNumber);
    end
end

if (verbose)
    toc;
end
%% Auxiliary functions
function fp=FractionalPart(x)
    fp=x-fix(x);



function rpara=ExMOGA(MOGAResults,index, plotf, splitn, ...
                      fitucf, fitdispf, fitbeta0f, fittunef, ...
                      plotdaf, tunesuc, betas0, tunes, saveOPAf, verbose, varargin)
% 
% Examines result of MOGA Scan, produces optics function plots and
% matches dispersion and achromat tunes 
%
% inputs:
%           MOGAResults: structure containig MOGA results (see MOGA.m
%           index     : index of line in "Values" in ScanResults structure from which 
%                       data will be retrieved. If negative the line number
%                       is used directly. Convenient to keep track of an
%                       individual regardless of he order of sorting o the
%                       ParetoFront array.
%           plotf     : if "Y", plots lattice functions for the extracted lattice
%                       as well for the descendent version with dispersion and tune fits
%           splitn    : breaks the elements in 2^splitn slices (if splitn is larger
%                       than 1)
%           fitucf    : if 'Y' fits the unit cell tunes to specified values
%           fitdispf  : if 'Y' matches dispersion to zero in long straight
%           fitbeta0f : if 'Y' fits the horizontal and vertical beta
%                       functions at centre of long stright
%           fittunef  : if "Y" fits lattice to match dispersion and fit period tunes
%                       all other fits.
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
%            verbose  : controls level of diagnostic printout
%
%            varargin{1}: DAOptions : alternative set of parameters for DA calcualtion.
%                        if not given , takes the values from the
%                        MOGAResults structure
% outputs:
%          rpara: structure containing lattice properties as produced by
%                 function atsummary_fast.m
%
%
%
LatticeOptData = MOGAResults.LatticeOptData;

isdipole    = LatticeOptData.isdipole;
nvars       = LatticeOptData.nvars;
LAT         = LatticeOptData.ACHRO;
UC          = LatticeOptData.UC;
Nper        = LatticeOptData.Nper;
Trb         = LatticeOptData.Trb;
scan_fams   = LatticeOptData.scan_fams;

RunNumber   = MOGAResults.RunNumber;
lb          = MOGAResults.lb;
ub          = MOGAResults.ub;

%
% check for backward compatibility
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

if(isfield(LatticeOptData,'lin_fams'))
    lin_fams = LatticeOptData.lin_fams;
else
    lin_fams = {};
end

%
% Parameters for dynamic aperture calculation
%% note that these may be different from the parameters at the time
%% MOGA was run and recorded in the structure 
%% MOGAResults.LatticeOptData.DAoptions
%
% If the argumetn DAOptions given to the function is an empty matrix, the
% corresponding DAOptions are taken from LatticeOptData structure in 
% MOGAResults. 
%
if (nargin>14)
    DAoptions     = varargin{1};
else
    DAoptions     = LatticeOptData.DAoptions;
end

chroms0      = DAoptions.chroms0; % Target chromaticity for one superperiod
TolChrom     = DAoptions.TolChrom;% Chromaticity tolerances
Nitchro      = DAoptions.Nitchro; % Max n. iterations of chromaticty correction
nturns       = DAoptions.nturns;  % number of turns
dp           = DAoptions.dp;      % dp/p
betax0       = DAoptions.betax0;  % hor beta for normalization
betay0       = DAoptions.betay0;  % ver beta for normalization

% DA calculation mode "border"
r0           = DAoptions.r0;      % initial guess [m];
nangs        = DAoptions.nang;    % number of angular steps
dang         = pi/(nangs-1);      % angular step size [rad}
res          = DAoptions.res;     % resolution [m]
DAalpha      = DAoptions.alpha;   % DA radius enlargement factor

%
% Default values for XmaxDA and YmaxDA are for use in scaling DA plots.
%
if(not(isfield(DAoptions,'XmaxDA')))
    XmaxDA     = 0.01;
else
    XmaxDA     = DAoptions.XmaxDA;  % Horizontal range is -Xmax to Xmax [m]
end

if(not(isfield(DAoptions,'YmaxDA')))
    YmaxDA     = 0.008;
else
    YmaxDA     = DAoptions.YmaxDA;  % Vertical range is 0 to Ymax [m]
end

if(not(isfield(DAoptions,'DAmode')))
    DAmode='border';
else
    DAmode = DAoptions.DAmode;
end

if (strcmp(DAmode,'grid'))
    % Parameters for "grid" dynamic aperture calculation
    
    npdax      = DAoptions.npdax;   % number of grid points in x direction is 2*npdax+1
    npday      = DAoptions.npday;   % number of grid points in y direction is  npday+1
    npDA       = DAoptions.npDA;    % total numbr of grid points
%
% Recalculates X0da and Y0da in case the data in DAoptions 
% does not come from the MOGAResults strucrure
%
    if (nargin>14)
        dx = XmaxDA/npdax; % grid stepsize in x [m]
        dy = YmaxDA/npday;   % grid stepsize in y [m]
        dxdy = dx*dy; % grid cell area [m**2]
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
    else
        dx         = DAoptions.dx;      % horizontal step size
        dy         = DAoptions.dy;      % vertical step size
        dxdy       = DAoptions.dxdy;    % grid cell area [m**2]
        X0da = DAoptions.X0da;    % horizontal coordinates of grid points [m]
        Y0da = DAoptions.Y0da;    % vertical coordinates of grid points [m]
    end
end

%
%
% Searches for desired individual in the Pareto Front
%
% Last column of ParetoFront table must contain index
%
indexcol = size(MOGAResults.ParetoFront,2); 
if (index<0)
    i=abs(index);
else
    try       
        i=find(MOGAResults.ParetoFront(:,indexcol)==index,1);
        if  isempty(i)
            if (strcmp(verbose,'Y')) 
                fprintf('Index %4d not found in ParetoFront; using i=1', index); 
            end
            i=1;
        end
    catch
        if (strcmp(verbose,'Y')) 
            fprintf('Index column not found in ParetoFront; using i=1');
        end
        i=1;
    end
end

Ks = MOGAResults.ParetoFront(i,1:nvars);

try
    iindex  = MOGAResults.ParetoFront(i,indexcol);
catch
    if (strcmp(verbose,'Y'))
        fprintf('index column not found in ParetoFront array \n');
    end
    iindex=1;
end

Knew = Ks;

if (strcmpi(verbose,'Y'))
    fprintf('**** \n');
    fprintf('%s Extracting individual n.%3d from MOGA set: %s \n', datetime, iindex, RunNumber);
    fprintf('**** \n');
    tic;
end

%
% Sets the decision variables for the retrieved individual
%

LAT  = setDVs(2, LAT,LatticeOptData, Ks);
UC   = setDVs(3, UC, LatticeOptData, Ks);

%loads linear lattice parameters for optimization mode 'SEXT'

if(strcmp(optMode,'SEXT'))
  DVLins=MOGAResults.DVLins;
  LAT = setLins(2,LAT,LatticeOptData,DVLins);
  UC  = setLins(3,UC,LatticeOptData,DVLins);
end
%
% Fits chromaticity for modes optMode='CHRO' and 'SEXT'
%
%
if(strcmp(optMode,'CHRO')||strcmp(optMode,'SEXT'))
     if (strcmp(verbose,'Y')) 
      fprintf('Fitting Chromaticity... \n'); 
     end
     try
        LAT_C=LAT; 
        k=0;    
        [~,chr]=tunechrom(LAT_C);
        while ( (abs(chr(1)-chroms0(1))>TolChrom(1))||...
                (abs(chr(2)-chroms0(2))>TolChrom(2)) &&(k<Nitchro))
           LAT_C=atfitchrom(LAT_C,chroms0,chrom_fams{1},chrom_fams{2});
           [~,chr]=tunechrom(LAT_C);
           k=k+1;
        end
        I_sd  = find(atgetcells(LAT_C, 'FamName', chrom_fams{1}));
        K_sd  = atgetfieldvalues(LAT_C, I_sd, 'PolynomB', {3});
        I_sfi = find(atgetcells(LAT_C, 'FamName', chrom_fams{2}));
        K_sfi = atgetfieldvalues(LAT_C, I_sfi, 'PolynomB', {3});
        SD  = K_sd(1);
        SFI = K_sfi(1);
     catch
        disp('Chromaticity fitting not possible');
        LAT_C = LAT;
     end
     LAT=LAT_C;
else
    SD=NaN;
    SFI=NaN;
end
%
% Plots lattice functions as derived from MOGA results directly
%
if (strcmp(plotf,'Y'))
    LAT_SP=LAT;
    if (splitn>1)
       for i=1:splitn
           LAT_SP=atinsertelems(LAT_SP,1:length(LAT_SP),0.5,[]);
        end
    end
    try
        lindata=atlinopt(LAT_SP,0.0,1:length(LAT_SP)+1);
        PlotBetaDisp(lindata,'MOGA bare');
    catch ME
        fprintf('Error in atlinopt for MOGA lattice \n');
        fprintf('Error message was:%s \n',ME.message);
    end
end

try
    rpara=atsummary_fast(LAT,isdipole);
catch ME
    fprintf('Error in atsummary full period from MOGA \n');
    fprintf('Error message was:%s \n',ME.message);
end
    
try
    rparaUC = atsummary(UC);
    tunesuc0 = rparaUC.tunes;
catch
    fprintf('error in unit cell tune calculation \n');
    tunesuc0 = [NaN NaN];
end

tunesuc1=tunesuc0;

LAT_UC=LAT;
Kold = Ks;

%
% Some linear lattice fits are only performed if mode is not "SEXT"
% and reverse bend mode is not U3 . To be improved to allow fits also in
% "U3" reverse bend mode.
%
if(not(strcmp(optMode,'SEXT')))
 if(strcmp(fitucf,'Y')&&not(isnan(tunesuc0(1)))&&not(isnan(tunesuc0(2)))&&(not(strcmp(revBmode,'U3'))))
    if (strcmp(verbose,'Y'))
        fprintf('Fitting unit cell tunes from %4.3f %4.3f to %4.3f %4.3f \n',...
             tunesuc0(1), tunesuc0(2), tunesuc(1), tunesuc(2) );
    end
    fam1 = 'reversebend_sim';
    if (strcmpi(optMode,'SIMP'))
        fam2 = 'sshdip';
    else
        fam2 = 'dip';
    end
%   UC_tune = atfittune(UC,  tunesuc, fam1, fam2);
%   UC_tune = atfittune(UC_tune,tunesuc , fam1, fam2);
    [UC_tune, its, penaltyuctune ftunes]= fittuneRS(UC, tunesuc, fam1, fam2, 1000, 1E-7, 'N');
%   [UC_tune, penaltyuctune] = fitfullTune(UC, fam1, fam2, lb(3:4), ub(3:4), tunesuc);
    try
        rpuc = atsummary(UC_tune);
        tunesuc1  = rpuc.tunes;
        if (strcmpi(verbose,'Y'))
            fprintf('Unit cell tunes fit to[ %8.5f %8.5f ]\n',...
                 tunesuc1(1), tunesuc1(2));
             fprintf (' in %5d iterations and penalty = %8.2e \n', ...
                 its, penaltyuctune );
        end
        Knew    = getDVs(3, UC_tune, LatticeOptData); 
        UC_tune = setDVs(3, UC_tune, LatticeOptData,Knew); 
        LAT_UC  = setDVs(2, LAT, LatticeOptData, Knew);
 %
 % Here Knew contains NaNs as  the unit cell does not have all families.
 % Below we fix this.
        Knew    = getDVs(2, LAT_UC, LatticeOptData);
                    
    catch ME
        fprintf('Error in atsummary of unit cell after unit cell fit \n');
        fprintf('Error message was:%s \n',ME.message);
        Dialog = questdlg('Revert lattice ?','Yes', 'No');
        switch Dialog
            case 'Yes'
                UC_tune=UC;
                LAT_UC = LAT;
                tunesuc1=tunesuc0;
                Knew = Kold;
                
            case 'No'
                
        end % switch
    end
 
    %
    % Checks lattice of the full period with the  updated unit cell
    %
    try
        rpara=atsummary_fast(LAT_UC,isdipole);
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
                fprintf('Lattice reverted \n');
            case 'No'
 
        end % switch    
    end
end

Kold=Knew;
LAT_md = LAT_UC;
   
if (strcmp(fitdispf,'Y'))
    if (strcmp(verbose,'Y'))
        fprintf('Matching dispersion \n');
    end
    switch optMode
        case {'SIMP','COMP'}
            rbfam = 'reversebendmc1_sim';
        case 'CHRO'
            switch revBmode
                case 'All'
                    rbfam = 'reversebendmc1_chro';
                case 'U3'
                    rbfam = 'qfm1_rbu3chro';
            end
    end
    Variab1 = atVariableBuilder(LAT_UC,...
               rbfam,...
               {'PolynomB',{1,2}},...
                lb(6),ub(6));
 
    Variab2 = atVariableBuilder(LAT_UC,...
               'reversebendmc2_sim',...
               {'PolynomB',{1,2}},...
               lb(7),ub(7));
           
%    Variables = [Variab1;Variab2];
    Variables = Variab1;
    Constr1  = atlinconstraint(1,{{'Dispersion',{1}}},0,0,1);
  
    [LAT_md, penalty_md, dmin_md] = atmatch(LAT_UC,Variables,Constr1,...
                         1E-8,1000,0,@fminsearch); %  @fmincon  
    
    if (strcmp(verbose,'Y'))
        fprintf('Dispersion matched with penalty = %6.2e \n', sum(penalty_md.^2));
    end

    % Makes sure K and PolynomB are the same to avoid warning messages
    % later
    
    Knew = getDVs(2,LAT_md,LatticeOptData);
    LAT_md = setDVs(2,LAT_md,LatticeOptData,Knew);
    
    
    try
        rpara=atsummary_fast(LAT_md,isdipole);
    catch ME
        fprintf('Error in atsummary of full period after dispersion fit \n');
        fprintf('Error message was:%s \n',ME.message);
        Dialog = questdlg('Revert lattice ?','Yes', 'No');
        switch Dialog
            case 'Yes'
                LAT_md = LAT_UC;
                Knew = Kold;
                fprintf('Lattice reverted \n');
                
            case 'No'

        end % switch 
    end
end

Kold = Knew;
LAT_bet = LAT_md;

if(strcmp(fitbeta0f,'Y'))
    if (strcmp(verbose,'Y'))
        fprintf('Fitting beta functions \n');
    end
    Variab1 = atVariableBuilder(LAT_bet,{'qfend_sim','qdend_sim'},...
                           {{'PolynomB',{1,2}},{'PolynomB',{1,2}}},...
                           {lb(1),lb(2)},{ub(1),ub(2)});
    Constr1 = atlinconstraint(1,{{'beta',{1}},{'beta',{2}}},...
                  betas0,betas0,[1,1]);
    
    [LAT_bet, penalty_bet, dmin_bet] = atmatch(LAT_md,Variab1,Constr1,...
                         1E-6,1000,0,@fminsearch); %  @fmincon  
    if (strcmp(verbose,'Y'))   
        fprintf('Betas matched with penalty = %6.2e \n', sum(penalty_bet.^2));
    end
    
    Knew = getDVs(2,LAT_bet,LatticeOptData);
    LAT_bet = setDVs(2,LAT_bet,LatticeOptData,Knew);
    
    try
        rpara=atsummary_fast(LAT_bet,isdipole);
    catch ME
        fprintf('Error in atsummary of full period after beta function fit \n');
        fprintf('Error message was: %s \n',ME.message);
        Dialog = questdlg('Revert lattice ?','Yes', 'No');
        switch Dialog
            case 'Yes'
                LAT_bet = LAT_md;
                Knew = Kold;
                fprintf('Lattice reverted \n');
                
            case 'No'
                
        end % switch 
    end    
end 
    
Kold = Knew;
LAT_tune = LAT_bet;
if(strcmp(fittunef,'Y'))
     qxfit = tunes(1)/Nper;
     qyfit = tunes(2)/Nper;
     if (strcmp(verbose,'Y'))
        fprintf('Fitting period tunes from [ %5.3f %5.3f ] to [ %5.3f %5.3f ] \n',...
              rpara.Qx_ring/Nper, rpara.Qy_ring/Nper,qxfit,qyfit );
     end
    
%     LAT_tune = atfittune(LAT_md,  [qxfit qyfit], 'qfend_sim', 'qdend_sim', 'UseIntegerPart');
%     LAT_tune = atfittune(LAT_tune,[qxfit qyfit] ,'qfend_sim', 'qdend_sim', 'UseIntegerPart');

     %lb = MOGAResults.lb(1:2);
     %ub = MOGAResults.ub(1:2);       
 %   [LAT_tune, penalty_tune] = fitfullTune(LAT_bet, 'qfend_sim',
 %   'qdend_sim', lb, ub, [qxfit qyfit]);
    switch optMode
        case {'COMP','SIMP'}
            tufam = {'qfend_sim';'qdend_sim'};
        case 'CHRO'
            tufam= {'qfend_rbu3chro';'qdend_rbu3chro'};
    end
     [LAT_tune, its, penalty_tune]= fittuneRS(LAT_bet, [qxfit qyfit],...
         tufam{1}, tufam{2}, 1000, 1E-7,'Y');
    if (strcmpi(verbose,'Y'))
        fprintf('Period Tune fit complete with penalty = %6.2e after %5d iterations \n', penalty_tune, its);
    end
 
    Knew = getDVs(2,LAT_tune,LatticeOptData);
    LAT_tune = setDVs(2,LAT_tune,LatticeOptData,Knew);
   
    try
        rpara=atsummary_fast(LAT_tune,isdipole);
        if (strcmp(verbose,'Y'))
            fprintf('Final ring tunes = [ %8.5f %8.5f ] \n', rpara.Qx_ring, rpara.Qy_ring);
        end
    
     catch ME
        fprintf('Error in atsummary of full period after tune fit \n');
        fprintf('Error message was: %s \n',ME.message);
        Dialog = questdlg('Revert lattice ?','Yes', 'No');
        switch Dialog
            case 'Yes'
                LAT_tune = LAT_bet;
                Knew = Kold;
                fprintf('Lattice reverted \n');
                
            case 'No'
                
        end % switch 
    end    
 end
else
    LAT_tune=LAT;
end

if (strcmp(plotf,'Y'))
    if (splitn>1)
        for i=1:splitn
            LAT_tune=atinsertelems(LAT_tune,1:length(LAT_tune),0.5,[]);
        end
    end
    try
        lindata=atlinopt(LAT_tune,0.0,1:length(LAT_tune)+1);
        PlotBetaDisp(lindata,'MOGA + fitted');
    catch ME
        fprintf('Error in atlinopt of final lattice \n');
        fprintf('Error message was: %s \n',ME.message);
    end
end

if (strcmp(plotdaf,'Y'))
      try
          if (strcmp(verbose,'Y'))   
                fprintf('Calculating Dynamic Aperture\n');
          end
          switch DAmode
              case 'border'
                    DAV = modelDA_sim_par(LAT_tune, r0, nangs, nturns, dp, 0.0, res, DAalpha);
                    figure;plot(DAV(:,1)*1000,DAV(:,2)*1000,'-ob');
                    xlabel('X [mm]'); ylabel('Y [mm]');grid;
                    xlim([-XmaxDA XmaxDA]*1000);ylim([0 YmaxDA]*1000);
                    rpara.DAV=DAV;
                    if (strcmpi(optMode,'CHRO')) 
                        DAV(:,1)=DAV(:,1)*sqrt(betax0/rpara.beta0(1));
                        DAV(:,2)=DAV(:,2)*sqrt(betay0/rpara.beta0(2));
                    end 
                    RDA = DAV.*DAV*1E6; % converts to mm**2
                    RA  = RDA(:,1)+RDA(:,2);
                    DA  = sum(RA)*dang/2;

              case 'grid'
                    DAV = calcDA_grid(LAT_tune, X0da, Y0da, nturns, dp);
                    DA  = sum(DAV)*dxdy;
                    if (strcmpi(LatticeOptData.optMode,'CHRO')) 
                        DA=DA*sqrt(betax0/rpara.beta0(1))*sqrt(betay0/rpara.beta0(2));
                    end 
                    DA=DA*1e6; % coverts to mm**2
                    DAM = zeros(npday+1,2*npdax+1);
                    k= 1;
                    for i=0:npday
                        for j= 1:2*npdax+1
                            DAM(npday+1-i,j)=DAV(k);
                            k=k+1;
                        end
                    end
                    DAM=DAM*255;
                    map=[0 0.75 0; 1 1 1];
                    figure;image([-XmaxDA*1000,XmaxDA*1000],[YmaxDA*1000,0],DAM);
                    ax=gca;
                    ax.YDir='normal';
                    colormap(map);
                    xlabel('X[mm]');
                    ylabel('Y[mm]');    
                    xlim([-XmaxDA XmaxDA]*1000);ylim([0 YmaxDA]*1000);grid;
            otherwise
                    fprintf('Uknown DA calculation mode %s \n', DAmode);
                    DA=NaN;
          end

      catch ME
        fprintf('Problems calculating dynamic aperture \n');
        fprintf('Error message was:%s \n',ME.message);
        DA=NaN;
      end
else
    DA=NaN;
end

% Adds info to output structure

rpara.Knew = Knew;
rpara.K0   = Ks;
rpara.tunesuc0=tunesuc0;
rpara.tunesuc1=tunesuc1;
rpara.index=iindex;
rpara.LatticeOptMode = LatticeOptData;
rpara.MOGARunNumber = RunNumber;
rpara.MOGAIndex=index;
rpara.SD=SD;
rpara.SFI=SFI;
rpara.DA=DA;

if(strcmp(optMode,'SEXT'))
    rpara.DVLins = DVLins;
end

%
% Saves OPA file
%
if(strcmp(saveOPAf,'Y'))
    switch revBmode
        case 'All'
          switch optMode
            case 'SIMP'
                filein = 'm4_Studies_InputfromMSOGA_AT_Template.opa';
                fileout = 'm4_Studies_InputfromMSOGA_AT.opa';
            case 'COMP'
                filein = 'm4_Studies_InputfromMSOGA_AT_Template_COMP.opa';
                fileout = 'm4_Studies_InputfromMSOGA_AT_COMP.opa';
            case 'CHRO'
                filein = 'm4_Studies_InputfromMSOGA_AT_Template_CHRO.opa';
                fileout = 'm4_Studies_InputfromMSOGA_AT_CHRO.opa';
            case 'SEXT'
                filein = 'm4_Studies_InputfromMSOGA_AT_Template_SEXT.opa';
                fileout = 'm4_Studies_InputfromMSOGA_AT_SEXT.opa';
            otherwise
                fprintf('Unknown optmization mode %s. Aborting \n', optMode);
            return
          end

        case 'U3'
           switch optMode
             case 'CHRO'
                filein = 'm4_Studies_InputfromMSOGA_AT_Template_CHROU3.opa';
                fileout = 'm4_Studies_InputfromMSOGA_AT_CHROU3.opa';

             otherwise
                   fprintf('Unknown optmization mode %s. Aborting \n', optMode);
                   return
           end
    end

    fileIDin = fopen(filein,'r');
    %pathnameout = 'C:\Users\pedtav\Documents\Accelerators\NewMachineStudies\Lattices\OrbitShift\OPA\';
    %fileIDout = fopen([pathnameout fileout],'w');
    fileIDout = fopen(fileout,'w');
    fprintf(fileIDout, '{ Decision Variables from Run  %s } \r\n', RunNumber);
    for i=1:nvars
       fprintf(fileIDout,'DV%02d = %8.5f ; { %s } \r\n',i,Knew(i),scan_fams{i});
    end
    %
    % if SEXT MODE includes linear optics parameters
    %
    if (strcmp(optMode,'SEXT'))
       nvars_lin=LatticeOptData.nvars_lin;
       for i=1:nvars_lin
           fprintf(fileIDout, 'DVLins%02d  = %8.5f ; { %s } \r\n',i,DVLins(i),lin_fams{i});
       end
    end
    %
    % in 'CHROM' or 'SEXT' modes, saves also the strengths of sd and sfi sexupoles, set
    % for chromaticty +1 in both planes
    %
    if (strcmp(optMode,'CHRO')||strcmp(optMode,'SEXT'))
       fprintf(fileIDout, 'SDMSOGA  = %8.5f ; \r\n', K_sd(1));
       fprintf(fileIDout, 'SFIMSOGA = %8.5f ; \r\n', K_sfi(1));
    end
    
    fprintf(fileIDout, 'RBKMSOGA = %8.5f ;\r\n', Trb);
    
    while 1
        tline = fgetl(fileIDin);
        if ~ischar(tline), break, end
        fprintf(fileIDout,'%s \r\n', tline);
     end
     fclose(fileIDin);
     fclose(fileIDout);
     fprintf('saved OPA file %s \n',fileout);
end

if (strcmp(verbose,'Y'))
    toc;
end

end

function fp=FractionalPart(x)
    fp=x-fix(x);
end



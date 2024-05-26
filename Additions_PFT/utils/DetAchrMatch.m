function DetMatchR = DetAchrMatch(varargin)
% DetAChrMatch
%   Returns AT Lattice cell array and beam parameters for ring composed of unit cells
%   with reverse bends and matching cells to dispersion-matched long
%   straights
%% inputs:
% manadatory inputs
%          LatticeOptData: structure containing the lattice optimization data -
% optional inputs 
%          fixrbk : if present reverse bending angle comes from
%                   LatticeOptData. Otherwise m4U is called with the rbk
%                   below as input
%          rbk     = Reverse bend angle  [rad]
%          verbose  : if present produces verbose output
%          plot     : if present produces lattice function plots
%          saveOPA  : if present saves OPA file with final strengths
%          nojoint  : if present, does not try to do joint dispersion and alpha
%                     matching
%          tunesuc  : [qx qy] target unit cell tunes
%          nittue   : number of iterations for tune mathcing
%          nitdisp  : number of iterations for dispersion matching
%          nitchro  : number of iterations for chromaticity matching
%          nitAlphas: number of iterations for alpha matching
%          TolX     : tolerance for changes in decision variables in
%                     when matching
%          TolTune  : tolerance for tune myachinh
%          TolDisp  : tolerance for dispersion matching
%          TolAlphas: tolearance for alpha matching
%          TolChro  : tolearance for chromaticty matching
%          lb: lower bounds for strengths to be used in fits; if empty uses
%              dx argument to calculate.
%          ub: lower bounds for strengths to be used in fits;  if empty uses
%              dx argument to calculate  
%          dx: relative (in percentage) range around initaal values for fits. 
%          chroms0 :  = [ksix, ksiy] target ring chromaticities
%
%% ouputs: structure with fields inputs and outputs 
%   dTheta = Reverse bend angle [mrad]
%   Krb    = Reverse Bend Quadrupole strength [m**-2]
%   Circ   = Circumference[m]
%   dCirC  = change in circumference compared to zero RB kick case [mm]
%   dfRF   = change in RF frequency [Hz]
%   TotAngle = total bend angle [deg]
%   yOffset  = vertical orbit offset at the position of the Reverse Bends [mm]
%   unpkick  = kick from Reverse bend due to the vertical offset above [mrad]
%   RbShift  = additional horizontal shift of Reverse bend in order to
%              achieve required total RB kick [mm]
%   emit     = emittance [pm rad]
%   Remit    = Ratio of emittance to unperturbed case [%]
%   Jx       = Horizontal damping partition
%   alphac   = momentum compaction [0.001]
%   sigma_E  = energy spread [0.001]
%   Qx, Qy   = Betatron Tunes
%   Ksix Ksiy = Chromaticities
%   betax0 betay0 disp0 = betatron and dispersion fucntions at center of
%   dipole [m]
%   Ksds  : strength of defocussign sextupole family [m**-3]
%   Ksfis : strength of focussing sextupooe family [m**-3]
%   RDTs  : Resonant Driving Terms (absolute values:
%   h11001,h00111),h20001,h00201,h10002,h10010,h10100,h21000,h30000,
%   h10110, h10020,h10200,h22000,h11110,h00220,h31000,h40000,h20110,
%   h11200,h20020,h20200,h00310,h00400,dnux_dJx,dnux_dJy, dnuy_dJy,
%   On-energy dynamic aperture:
%   apx = positive horizontal [mm]
%   amx = negative horizontal [mm]
%   apy = vertical
%   acpx = apx^2/betax0 [mm]
%   acmx = amx^2/betax0 [mm]
%   acy  = apy^2/betay0 [mm] 
%   ac   = (acpx+acmx)*acy [mm**2]
%   qx   = horizontal cell tune
%   qy   = vertical cell tune
%   nuc  = number of unit cells
%% Input argument parsing
LatticeOptData = getargs(varargin,[]);
plotf          = any(strcmpi(varargin,'plot'));
fixedThetaf    = any(strcmpi(varargin,'fixrbk'));
fittunef       = any(strcmpi(varargin,'fittune'));
fitchromf      = any(strcmpi(varargin,'fitchrom')); 
plotdaf        = any(strcmpi(varargin,'plotda')); 
savef          = any(strcmpi(varargin,'save'));

saveOPAf       = any(strcmpi(varargin,'saveOPA'));
verbosef       = any(strcmpi(varargin,'verbose'));
nojointf       = any(strcmpi(varargin,'nojoint'));

dx             = getoption(varargin,'dx',10);
lb             = getoption(varargin,'lb',[]);
ub             = getoption(varargin,'ub',[]);
dTheta         = getoption(varargin,'rbk',0.0);
verboselevel   = getoption(varargin,'verboselevel',0);
TolX           = getoption(varargin,'TolX',1E-8);

tunes          = getoption(varargin,'tunes',[46.20 16.28]);
tunesuc        = getoption(varargin,'tunesuc',[3/7 1/7]);
nittune        = getoption(varargin,'nittune',LatticeOptData.nittune);
nitdisp        = getoption(varargin,'nitdisp',LatticeOptData.nitdisp);
nitAlpha       = getoption(varargin,'nitAlpha',LatticeOptData.nitAlpha);
TolTune        = getoption(varargin,'TolTune',LatticeOptData.TolTune);
TolDisp        = getoption(varargin,'TolDisp',LatticeOptData.TolDisp);
TolAlpha       = getoption(varargin,'TolAlpha',LatticeOptData.TolAlpha);
DAoptions      = getoption(varargin,'DAoptions',LatticeOptData.DAoptions);

chroms0        = getoption(varargin,'chroms0',DAoptions.chroms0);
Nitchro        = getoption(varargin,'nitchro',DAoptions.Nitchro);
TolChrom       = getoption(varargin,'TolChrom',DAoptions.TolChrom);

%% preamble
if (verbosef)
    fprintf('******* \n');
    fprintf('%s Determinist Achromat match for unit cell tune = %5.3f %5.3f \n', datetime, tunesuc);
    fprintf('******* \n');
end
UC            = LatticeOptData.UC;
IMC1          = LatticeOptData.IMC1;
HACHRO        = LatticeOptData.HACHRO;
ACHRO         = LatticeOptData.ACHRO;
ACHRO_fit     = ACHRO;
IfamsF        = LatticeOptData.IfamsF;
IfamsUC       = LatticeOptData.IfamsUC;
IfamsIMC1     = LatticeOptData.IfamsIMC1;
I_famsAllF    = LatticeOptData.IfamsAllF;
I_famsAllUC   = LatticeOptData.IfamsAllUC;
nvars         = LatticeOptData.nvars;
lattMode      = LatticeOptData.lattMode;
optMode       = LatticeOptData.optMode;
Nper          = LatticeOptData.Nper;

scan_fams     = LatticeOptData.scan_fams;
uctune_fams   = LatticeOptData.uctune_fams;
disp_fams     = LatticeOptData.disp_fams;
betaxy0_fams  = LatticeOptData.betaxy0_fams;
alpha0_fams   = LatticeOptData.alpha0_fams;
ringtune_fams = LatticeOptData.ringtune_fams;
chrom_fams    = LatticeOptData.chrom_fams;
isdipole      = LatticeOptData.isdipole;

if (isfield(LatticeOptData,'lin_fams'))
    lin_fams      = LatticeOptData.lin_fams;
else
    lin_fams={};
end
All_fams      = LatticeOptData.All_fams;
nallfams      = LatticeOptData.nallfams;
stdfamlist    = LatticeOptData.All_stdfamlist; 
nstdfamlist   = LatticeOptData.All_nstdfamlist;

% Iuctune = zeros(size(uctune_fams,1));
% for i=1:size(uctune_fams,1)
%     Iuctune(i)=find(strcmp(All_fams,uctune_fams{i}));
% end
% Idisp = zeros(size(disp_fams,1));
% for i=1:size(size(disp_fams,1))
%     Idisp(i)=find(strcmp(All_fams,disp_fams{i}));
% end
% Ibetaxy0 = zeros(size(betaxy0_fams,1));
% for i=1:size(betaxy0_fams,1)
%     Ibetaxy0(i) = find(strcmp(All_fams,betaxy0_fams{i}));
% end
% Iringtune = zeros(size(ringtune_fams,1));
% for i=1:size(ringtune_fams,1)
%     Iringtune(i)=find(strcmp(All_fams,ringtune_fams{i}));
% end
% Idvs = zeros(LatticeOptData.nvars);
% for i=1:nvars
%     Idvs(i)=find(strcmp(All_fams,scan_fams{i}));
% end
% Ilin = zeros(nlinfams);
% for i=1:nlinfams
%      Ilin(i)=find(strcmp(All_fams,lin_fams{i}));
% end

X0 = getAllfams(2,ACHRO,LatticeOptData);

if (not(size(lb,2)==nallfams)||not(size(ub,2)==nallfams))
    lb = X0-abs(X0)*dx/100;
    ub = X0+abs(X0)*dx/100;
end

if (verbosef)
    fprintf('lower bounds set to: ');
    fprintf('%5.3f ', lb(1:nallfams));
    fprintf('\r\n');
    fprintf('upper bounds set to: ');
    fprintf('%5.3f ', ub(1:nallfams));
    fprintf('\r\n');
end

% Kall = zeros(LatticeOptData.nallfams,1);
% for i=1:nvars
%     Kall(Idvs(i))=Ks(i);
% end
% for i=1:nlinfams
%     Kall(Ilin(i))=MOGAResults.DVLins(i);
% end


frf=99.931E6;
RbLength = 0.15;

if (fixedThetaf)
    dTheta = LatticeOptData.Trb;
else
    fprintf('Calling m4U to redefine reverse bend angle to %5.2f mrad\n', dTheta*1000);
    m4U(lattMode,optMode,dTheta);
end   

%% Orbit geometry
if (verbosef)
    disp('Calculating orbit geometry...');
    tic;
end

[x2d y2d a2d] = Survey2D(UC,1.5*pi/180);
yOffset = 6.509501134666746-y2d(round(length(y2d)/2))*1000;
if (verbosef)
    toc;
end
Penalty = 0.0;
%% Matches unit cell tunes
%
try
    ringpara = atsummary(UC);
    tunesuc0 = ringpara.tunes;
catch ME
    fprintf('Error in DetAchrMatc:unit cell tune calculation \n');
    fprintf('Error message was:%s \n',ME.message);
    ringpara=[];
    return
end
tunesuc1=tunesuc0;
X0_old = X0;
X0_new = X0_old;
%
UC_tune= UC;
if(not(isnan(tunesuc0(1)))&&not(isnan(tunesuc0(2)))&&not(isnan(tunesuc(1)))&&not(isnan(tunesuc(2))))
    if (verbosef)
        fprintf('Fitting unit cell tunes from %4.3f %4.3f to %4.3f %4.3f \n',...
             tunesuc0(1), tunesuc0(2), tunesuc(1), tunesuc(2) );
        tic;
    end
    
    try
        [UC_tune, its, penaltyuctune, ftunes]= fittuneRS(UC, tunesuc, uctune_fams{1}, uctune_fams{2},'maxits',nittune,'Tol',TolTune,'UseIntegerPart',false);
        Penalty = sqrt(Penalty^2 + penaltyuctune^2);
        if (verbosef)
            fprintf('Unit cell tunes fit to[ %8.5f %8.5f ]\n',...
                 tunesuc1(1), tunesuc1(2));
             fprintf (' in %5d iterations and penalty = %8.2e \n', ...
                 its, penaltyuctune );
        end
        X0_new    = getAllfams(3,UC_tune,LatticeOptData);
        UC_tune   = setAllfams(3,UC_tune,LatticeOptData,X0_new); % to guarantee PolynomB and K agree
        ACHRO_fit = setAllfams(2,ACHRO_fit,LatticeOptData,X0_new);
        X0_new    = getAllfams(2,ACHRO_fit,LatticeOptData);
        IMC1      = setAllfams(4,IMC1,LatticeOptData,X0_new);
        X0_old    = X0_new;
    catch ME
        fprintf('Error in DetAchrMatc:unit cell tune natching \n');
        fprintf('Error message was:%s \n',ME.message);
        ringpara=[]; 
    end
    try
        ringpara  = atsummary(UC_tune);
        tunesuc1  = ringpara.tunes;
    catch ME
        fprintf('Error in DetAchrMatch:atsummary of unit cell after unit cell fit \n');
        fprintf('Error message was:%s \n',ME.message);
        Dialog = questdlg('Revert lattice ?','Yes', 'No');
        switch Dialog
            case 'Yes'
                UC_tune=UC;
                ACHRO_fit = ACHRO;
                tunesuc1=tunesuc0;
                X0_new = X0_old;
            case 'No'
        end % switch
    end
    if (verbosef)
        toc;
    end
    if (plotf)
       try
          figure;atplot(UC_tune);title('Tune Matched')
       catch ME
          fprintf('Error in atplot of tune fitted Unit cell  \n');
          fprintf('Error message was: %s \n',ME.message);
          return
       end
    end
end

%% Propagates optics into the dispersion suppressor/matching cell
%
[twiss_uc,~,~] = atlinopt(UC_tune,0,1);
%
%% Matches dispersion
%
if (verbosef)
    disp('Matching Dispersion ...');
    tic;
end 

for i=1:length(LatticeOptData.DF_stdfamlist)
    kf = find(strcmp(All_fams,disp_fams{LatticeOptData.DF_stdfamlist(i)}));
    Variables(i) = atVariableBuilder(IMC1,disp_fams{LatticeOptData.DF_stdfamlist(i)},...
                   {'PolynomB',{1,2}},X0_new(kf),lb(kf),ub(kf)); 
end

for i=1:length(LatticeOptData.DF_nstdfamlist)
    kf = find(strcmp(All_fams,disp_fams{LatticeOptData.DF_nstdfamlist(i)}));
    Variables(i+length(LatticeOptData.DF_stdfamlist)) = ...
                 atVariableBuilder(@(R,K)setAllfams(4,R,LatticeOptData,...
                 cat(2,nan(1,kf-1),K,nan(1,nallfams-kf))),...
                 X0_new(kf),lb(kf),ub(kf));                         
end

%Variab1a = atVariableBuilder(IMC1,...
%                             disp_fams{1},...
%                             {'PolynomB',{1,2}});
%Variab1b = atVariableBuilder(IMC1,...
%                             disp_fams{2},...
%                             {'PolynomB',{1,2}});

%Variab2 = atVariableBuilder(@(R,K)setLinMatch(2,R,LatticeOptData,...
%                            cat(2,nan(1,nstdfamlist(i)-1),K,nan(1,nLMfams-nstdfamlist(i)))),...
%                            X0(nstdfamlist(i)),lb(nstdfamlist(i)),ub(nstdfamlist(i)));

%Variab2 = atVariableBuilder(@(R,K)setAllfams(4,R,LatticeOptData,...
%                            [NaN NaN NaN NaN K NaN NaN NaN NaN NaN]),...
%                            X0_new(5),lb(5),ub(5));

%Variab2 = atVariableBuilder(@(R,K)setDVs(4,R,LatticeOptData,...
%                            [NaN NaN NaN NaN K NaN NaN]),...
%                            DVs(4),lb(4),ub(4));

%Variables = [Variab1a, Variab1b, Variab2];
%Variables = [Variab1a, Variab1b];

Constraints  = [];
Constraints  = atlinconstraint(length(IMC1)+1,...
                          {{'Dispersion',{1}},...
                           {'Dispersion',{2}}},...
                          [0 0], [0 0], [1 1]);

%Constraints  = atlinconstraint(length(IMC1)+1,...
%                          {{'Dispersion',{1}}},...
%                          [0], [0], [1]);

try 
    [IMC1_md, penalty_md, dmin_md] = atmatch_mod(IMC1,Variables,Constraints,...
                         TolDisp,TolX,nitdisp,verboselevel,@fminsearch, ...
                         twiss_uc);
                     
    if (verbosef)
        toc;
        fprintf('Dispersion matched with penalty = %8.3e \n', sum(penalty_md.^2));
    end
    Penalty    = sqrt(Penalty^2  + sum(penalty_md.^2));
    X0_new     = getAllfams(4,IMC1_md,LatticeOptData);
    IMC1_md    = setAllfams(4,IMC1_md,LatticeOptData,X0_new); % to guarantee PolynomB and K agree
    ACHRO_fit  = setAllfams(2,ACHRO_fit,LatticeOptData,X0_new);
    X0_new     = getAllfams(2,ACHRO_fit,LatticeOptData);
                 
catch ME
        fprintf('Error in DetAchrMatch:dispersion matching. Aborting... \n');
        fprintf('Error message was:%s \n',ME.message);
        return
end

if (plotf)
    try
      figure;atplot(IMC1_md,'twiss_in',twiss_uc);title('Dispersion Matched')
    catch ME
        fprintf('Error in atplot of Dispersion fitted lattice \n');
        fprintf('Error message was: %s \n',ME.message);
        return
    end
end

%% Matches betatron function derivatives
%
if (verbosef)
    fprintf('Matching alphas...\n');
    tic;
end

for i=1:length(LatticeOptData.ALP_stdfamlist)
    kf = find(strcmp(All_fams,alpha0_fams{LatticeOptData.ALP_stdfamlist(i)}));
    Variables(i) = atVariableBuilder(IMC1,alpha0_fams{LatticeOptData.ALP_stdfamlist(i)},...
                   {'PolynomB',{1,2}},X0_new(kf),lb(kf),ub(kf)); 
end

for i=1:length(LatticeOptData.ALP_nstdfamlist)
    kf = find(strcmp(All_fams,alpha0_fams{LatticeOptData.ALP_nstdfamlist(i)}));
    Variables(i+length(LatticeOptData.ALP_stdfamlist)) = ...
                 atVariableBuilder(@(R,K)setAllfams(4,R,LatticeOptData,...
                 cat(2,nan(1,kf-1),K,nan(1,nallfams-kf))),...
                 X0_new(kf),lb(kf),ub(kf));                         
end

%Variab3 = atVariableBuilder(IMC1_md,...
%                            {scan_fams{1},scan_fams{2}},...
%                            {{'PolynomB',{1,2}},{'PolynomB',{1,2}}},...
%                            {lb(1),lb(2)},{ub(1),ub(2)});
Constraints = atlinconstraint(length(IMC1_md)+1,...
                          {{'alpha',{1}},...
                           {'alpha',{2}}},...
                          [0 0], [0 0], [1 1]);
try
    [IMC1_mdb,penalty_mdb,dmin_mdb] = atmatch_mod(IMC1_md,Variables,Constraints,...
                     TolAlpha,TolX,nitAlpha,verboselevel,@fminsearch,twiss_uc);
    if (verbosef)
        fprintf('Alphas matched with penalty = %8.3e \n',sum(penalty_mdb.^2));
    end
    Penalty    = sqrt(Penalty^2  + sum(penalty_mdb.^2));
    X0_new     = getAllfams(4,IMC1_mdb,LatticeOptData);
    IMC1_mdb   = setAllfams(4,IMC1_mdb,LatticeOptData,X0_new); % to guarantee PolynomB and K agree
    ACHRO_fit  = setAllfams(2,ACHRO_fit,LatticeOptData,X0_new);
    X0_new     = getAllfams(2,ACHRO_fit,LatticeOptData);
catch ME
    fprintf('Error in DetAchrMatch:dispersion matching. Aborting... \n');
    fprintf('Error message was:%s \n',ME.message);
    return
end

if (plotf)
    try
      figure;atplot(IMC1_mdb,'twiss_in',twiss_uc);title('Alphas Matched')
    catch ME
        fprintf('Error in atplot of Alphas fitted lattice \n');
        fprintf('Error message was: %s \n',ME.message);
        return
    end
end

%
% if fit is not accepted, then try a joint fit of dispersion and alphas
%
%IMC1_mdbj = IMC1_mdb;

% only checks for joint fits if nojointf is 'N'
%if (not(strcmp(nojointf,'Y')))
%    disjointfit = questdlg('Keep this fit ?');
%else
%    disjointfit = 'Yes';
%end

%if(not(strcmp(disjointfit,'Yes')))
%
%% Alternative match of dispersion and alphas simultaneously
%
%    if (strcmp(verbose,'Y'))
%        tic;
%        fprintf('Joint matching of dispersion and alphas \n');
%    end
% 
%     Variables   = [Variab1a,Variab1b,Variab2,Variab3];
%     Constraints = [Constr1,Constr2];
% 
%     [IMC1_mdbj, penalty_mdbj, dmin_mdbj] = atmatch_mod(IMC1_T,Variables,...
%                          Constraints,TolX, TolDisp,nit,...
%                          verboselevel,@fminsearch,...
%                          twiss_uc);
%     if (strcmp(verbose,'Y'))
%         toc;
%         fprintf('Alphas and dispersion matched with penalty = %8.3e \n', sum(penalty_mdbj.^2));
%         disp(dmin_mdbj);
%     end
%     if (strcmp(plotf,'Y'))
%         lindata=atlinopt(IMC1_mdbj,0.0,1:length(IMC1_mdbj)+1,... 
%         'twiss_in', twiss_uc);
%         PlotBetaDisp(lindata,'Dispersion+Alpha Joint Match');
%     end
% end

%twiss_0 = atlinopt(IMC1_mdbj, 0.0, length(IMC1_mdbj) + 1, 'twiss_in', twiss_uc);
%MC1 = flip(IMC1_mdbj);

%
% sets the new quad strengths to the whole achromat
%
% DVsMC  = getDVs(4,IMC1_mdb,LatticeOptData);
% DVsUC  = getDVs(3, UC_T, LatticeOptData);
% ACHRO  = setDVs(2,ACHRO,LatticeOptData,DVsMC);
% ACHRO  = setDVs(2,ACHRO,LatticeOptData,DVsUC);
% DVs    = getDVs(2,ACHRO,LatticeOptData);

%
%Zeros all sextupoles
%
%ACHRO  = scalesext(ACHRO,'sd',0);
%ACHRO  = scalesext(ACHRO,'sfm',0);
%ACHRO  = scalesext(ACHRO,'sfi',0);
%ACHRO  = scalesext(ACHRO,'sfo',0);
%ACHRO  = scalesext(ACHRO,'sdend',0);

%RINGRB   = repmat(ACHRORB,20,1);

%lindata_mc=atlinopt(MC1RB,0.0, 1:length(MC1RB)+1, 'twiss_in', twiss_0);

%lindata_mc(length(MC1RB)+1)
%PlotBetaDisp(lindata_mc);
%lindata_achrb=atlinopt(ACHRORB,0.0, 1:length(ACHRORB)+1, 'twiss_in', twiss_0);
%PlotBetaDisp(lindata_achrb);

if (plotf)
    try
      figure;atplot(ACHRO_fit);title('Full Achromat');
    catch ME
        fprintf('Error in atplot of Full achromat lattice \n');
        fprintf('Error message was: %s \n',ME.message);
        return
    end
end

%% Matches whole ring tunes
ACHRO_old = ACHRO_fit;
if(fittunef)
     qxfit = tunes(1)/Nper;
     qyfit = tunes(2)/Nper;
     if (verbosef)
        fprintf('Fitting period tunes from [ %5.3f %5.3f ] to [ %5.3f %5.3f ] \n',0.0,0.0,qxfit,qyfit);
%              rpara.Qx_ring/Nper, rpara.Qy_ring/Nper,qxfit,qyfit )
                      
     end
    
     [ACHRO_fit, its, penalty_tune]= fittuneRS(ACHRO_fit, [qxfit qyfit],...
         ringtune_fams{1}, ringtune_fams{2}, 'maxits', nittune, 'Tol', TolTune,'UseIntegerPart', true);
    if (verbosef)
        fprintf('Period Tune fit complete with penalty = %6.2e after %5d iterations \n', penalty_tune, its);
    end
 
    Penalty    = sqrt(Penalty^2  + sum(penalty_tune.^2));
    X0_new     = getAllfams(2,ACHRO_fit,LatticeOptData);
    ACHRO_fit  = setAllfams(2,ACHRO_fit,LatticeOptData,X0_new);
    X0_new     = getAllfams(2,ACHRO_fit,LatticeOptData);

    try
        rpara=atsummary_fast(ACHRO_fit,isdipole);
        if (verbosef)
            fprintf('Final ring tunes = [ %8.5f %8.5f ] \n', rpara.Qx_ring, rpara.Qy_ring);
        end
        if (plotf)
          figure;atplot(ACHRO_fit);title('Fnal Linear Match')
        end

     catch ME
        fprintf('Error in atsummary of full period after tune fit \n');
        fprintf('Error message was: %s \n',ME.message);
        Dialog = questdlg('Revert lattice ?','Yes', 'No');
        switch Dialog
            case 'Yes'
                ACHRO_fit = ACHRO_old;
                X0_new=X0_old;
                fprintf('Lattice reverted \n');
                
            case 'No'

        end % switch 
    end    
end
%% Matches chromaticity 
%
if(fitchromf&&~isempty(chrom_fams))
    ACHRO_old=ACHRO_fit;
     if (verbosef) 
      fprintf('Fitting Chromaticity... \n'); 
     end
     try 
        [ACHRO_fit, Penalty_chro, its]=fitchroit(ACHRO_fit, chrom_fams, chroms0, Nitchro, TolChrom); 
        I_sc1  = find(atgetcells(ACHRO_fit, 'FamName', chrom_fams{1}));
        K_sc1  = atgetfieldvalues(ACHRO_fit, I_sc1, 'PolynomB', {3});
        I_sc2  = find(atgetcells(ACHRO_fit, 'FamName', chrom_fams{2}));
        K_sc2  = atgetfieldvalues(ACHRO_fit, I_sc2, 'PolynomB', {3});
        Sc1    = K_sc1(1);
        Sc2    = K_sc2(1);
        Penalty= sqrt(Penalty^2  + sum(Penalty_chro.^2));
        if (verbosef)
            fprintf('Chromaticity matched with penalty = %6.2e in %2d iterations\n', Penalty, its);
        end
            X0_new    = getAllfams(2,ACHRO_fit,LatticeOptData);
            ACHRO_fit = setAllfams(2,ACHRO_fit,LatticeOptData,X0_new);
     catch ME
        fprintf('Error in DetAChroMatch: chromaticity fit \n');
        fprintf('Error message was: %s \n',ME.message);
        ACHRO_fit = ACHRO_old;
     end
     
     try
        rpara=atsummary_fast(ACHRO_fit,isdipole);
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
    DA=CalcPlotDA(ACHRO_fit,DAoptions,'plot');
    fprintf('Dynamic Aperture = %4.2f mm**2 \n',DA);
else
    DA=NaN;
end
%% Collects data
if (verbosef)
    disp('Calculating optics parameters for the final lattice...');
    tic;
end
try
    ringpara=atsummary_fast(ACHRO_fit,LatticeOptData.isdipole);
catch ME
    fprintf('Error in atsummary of Full achromat lattice \n');
    fprintf('Error message was: %s \n',ME.message);
    ringpara.circ=NaN;
    ringpara.naturalEmittance=NaN;
    ringpara.naturalEnergySpread=NaN;
    ringpara.tunes=[NaN NaN];
    ringpara.damping = [NaN NaN NaN];
    ringpara.chromaticity= [NaN NaN];
    ringpara.alphac = NaN;
    ringpara.Qx_ring = NaN;
    ringpara.Qy_ring = NaN;
    ringpara.beta0 = [NaN NaN];
    ringpara.etax0 = NaN;
end
if(verbosef)
    toc;
end
dCirc = (ringpara.circ-528.0)*1000;
dfRF  = dCirc/528.0*frf/1000;
%unpkick = RbLength*Krb*yOffset;
%RbShift = (dTheta*1000 - unpkick)/RbLength/Krb;

Ksixring = ringpara.chromaticity(1)*20;
Ksiyring = ringpara.chromaticity(2)*20;

DetMatchR.inputs.Trb             = dTheta;
DetMatchR.inputs.nojointfit      = nojointf;
DetMatchR.inputs.LatticeOptData  = LatticeOptData;
DetMatchR.inputs.tunesuc         = tunesuc;
DetMatchR.inputs.tunes           = tunes;
DetMatchR.inputs.TolX            = TolX;
DetMatchR.inputs.nittune         = nittune;
DetMatchR.inputs.nitAlpha        = nitAlpha;
DetMatchR.inputs.TolTune         = 'TolTune';
DetMatchR.inputs.TolDisp         = TolDisp;
DetMatchR.inputs.TolAlpha        = TolAlpha;

DetMatchR.inputs.lb              = lb;
DetMatchR.inputs.ub              = ub;
DetMatchR.inputs.fittunef        = fittunef;
DetMatchR.inputs.fitchromf       = fitchromf; 
DetMatchR.inputs.nojointf        = nojointf;
DetMatchR.inputs.dx              = dx;
DetMatchR.inputs.lb              = lb;
DetMatchR.inputs.ub              = ub;
DetMatchR.inputs.tunes           = tunes; 

DetMatchR.inputs.chroms0         = chroms0;

DetMatchR.inputs.Nitchro        = Nitchro;
DetMatchR.inputs.TolChrom       = TolChrom;

DetMatchR.outputs.ringpars       = ringpara;
DetMatchR.outputs.ringpars.dfRF  = dfRF;
DetMatchR.outputs.ringpars.yOffset  = yOffset;
%DetMatchR.outputs.ringpara.unpkick  = unpkick;
%DetMatchR.outputs.ringpara.RbShitf  = RbShift;
DetMatchR.outputs.ringpars.Ksixring = Ksixring;
DetMatchR.outputs.ringpars.Ksiyring = Ksiyring;
DetMatchR.outputs.Penalty = Penalty;
DetMatchR.outputs.XAll  = X0_new;
DetMatchR.outputs.ACHROMAT= ACHRO_fit;


%% Saves OPA file
if(saveOPAf)
   RunNumber = strcat('DeterministicMatch_',datestr(now,30));
   saveOPA(X0_new,[],[],LatticeOptData,RunNumber);
end
%% Saves results file
if(savef)
    filename = strcat('DetMatch_',datestr(now,30));
    DetMatchR.RunNumber = filename;
    fprintf('Saving file %s \n', filename);
    filename=strcat('DeterministicMatch/',filename);
    try
        save(filename,'DetMatchR');
    catch
        fprintf('%s Problems saving Deterministic Match Results file', datetime)
    end
end
          


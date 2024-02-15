function rp = RBAchroTuneParam_AT2(LatticeOptData, fixedThetaf, ...
                 dTheta, verbosef, verboselevel, plotf, saveOPAf, ...
                 nojointf, uctunes, nit, TolTunes, TolDisp, ...
                 TolAlphas, TolX, lb, ub, varargin)
% RBParam 
%   Returns beam parameters for ring composed of unit cells
%   with reverse bends and matching cells to dispersion-matched long
%   straights
%   Input:
%          LatticeOptData: structure containg the lattice optization data -
%          fixedThetaf : if 'Y' reverse bending angle comes from
%                        LatticeOptData
%          dTheta = Reverse bend angle  [rad]
%          verbosef : if 'Y' produces verbose output
%          plotf    : if yes produces lattice function plots
%          saveOPAf : saves OPA file with final strengths
%          nojointf : if 'Y', does not try to do joint dispersion and alpha
%                     matching
%          uctunes  : [qx qy] target unit cell tunes
%          nit      : number of iterations for fits
%          TolX     : tolerance on decision variables
%          TolDisp  : tolerance for dispersion matching
%          TolAlphas: tolearance for alpha matching
%          lb: lower bounds for strengths to be used in fits (1 to 7)
%          ub: lower bounds for strengths to be used in fits (1 to 7)
%          varargin{1} = targetchrom = [ksix, ksiy] target ring chromaticties
%
%   Ouput: structure with fields inputs and outputs 
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

UC        = LatticeOptData.UC;
IMC1      = LatticeOptData.IMC1;
HACHRO    = LatticeOptData.HACHRO;
ACHRO     = LatticeOptData.ACHRO;
scan_fams = LatticeOptData.scan_fams;
IfamsF    = LatticeOptData.IfamsF;
IfamsUC   = LatticeOptData.IfamsUC;
IfamsIMC1 = LatticeOptData.IfamsIMC1;
nvars     = LatticeOptData.nvars;

frf=99.931E6;
RbLength = 0.15;

if (nargin>16)
    chroms = varargin{1};
    fitchrom='Yes';
else
    fitchrom='No';
end

if (strcmp(fixedThetaf,'Y'))
    dTheta = LatticeOptData.Trb;
else
    max4_UpgradeStudies_20230305_AT2(1,'COMP','Yes',dTheta)
end   

if (strcmp(verbosef,'Y'))
    disp('Calculating orbit geometry...');
    tic;
end

[x2d y2d a2d] = Survey2D(UC,1.5*pi/180);
yOffset = 6.509501134666746-y2d(round(length(y2d)/2))*1000;
if (strcmp(verbosef,'Y'))
    toc;
end
%
%fits desired tunes
%
if (strcmp(verbosef,'Y'))
    disp('Fitting Tunes...');
    tic;
end
UC_T = UC;
%UCRB_T=atfittune(UCRB,tunes,'reversebend','dip');
%UCRB_T=atfittune(UCRB_T,tunes,'reversebend','dip');
try
  [UC_T, its, penalty, ftunes]=fittuneRS(UC,uctunes,scan_fams{3},scan_fams{4},...
                                 nit, TolTunes,'No');
catch
    fprintf('Problems fitting tunes %8.3f %8.3f \n', uctunes(1),uctunes(2));
end
%
% finds the new quadrupole strengths in the unit cell
%
DVs = getDVs(3,UC_T,LatticeOptData);
Krb = DVs(3);

if (strcmp(verbosef,'Y'))
    fprintf('Tunes fitted to %10.5f %10.5f with \n', ftunes(1),ftunes(2));
    fprintf('Revbend = %8.3f \n', DVs(3)); 
    fprintf('Dipole  = %8.3f \n', DVs(4));
end

%
% sets updated quadrupole strengths in the matching cell and achromat
%
IMC1_T = setDVs(4,IMC1,LatticeOptData,DVs);

% 
% propagates optics into the dispersion suppressor/matching cell
%
[twiss_uc,~,~] = atlinopt(UC_T,0,1);
%
% Matches dispersion
%
if (strcmp(verbosef,'Y'))
    disp('Matching Dispersion ...');
    tic;
end

Variab1a = atVariableBuilder(IMC1,...
                             scan_fams{6},...
                             {'PolynomB',{1,2}},...
                             lb(6),ub(6));
Variab1b = atVariableBuilder(IMC1,...
                             scan_fams{7},...
                             {'PolynomB',{1,2}},...
                             lb(7),ub(7));

Variab2 = atVariableBuilder(@(R,K)setDVs(4,R,LatticeOptData,...
                            [NaN NaN NaN NaN K NaN NaN]),...
                            DVs(4),lb(4),ub(4));

%Variables = [Variab1a, Variab1b, Variab2];
Variables = [Variab1a, Variab1b];

Constr1  = atlinconstraint(length(IMC1)+1,...
                          {{'Dispersion',{1}},...
                           {'Dispersion',{2}}},...
                          [0 0], [0 0], [1 1]);

[IMC1_md, penalty_md, dmin_md] = atmatch_mod(IMC1_T,Variables,Constr1,...
                         TolX, TolDisp,100,verboselevel,@fminsearch, ...
                         twiss_uc);
                     
if (strcmp(verbosef,'Y'))
    toc;
    fprintf('Dispersion matched with penalty = %8.3e \n', sum(penalty_md.^2));
    disp(dmin_md);
end

if (strcmp(plotf,'Y'))
    lindata=atlinopt(IMC1_md,0.0,1:length(IMC1_md)+1, 'twiss_in', twiss_uc);
    PlotBetaDisp(lindata,'Dispersion Match');
end
%
% Matches betatron function derivatives
%
if (strcmp(verbosef,'Y'))
    disp('Matching alphas...');
    tic;
end

Variab3 = atVariableBuilder(IMC1_md,...
                            {scan_fams{1},scan_fams{2}},...
                            {{'PolynomB',{1,2}},{'PolynomB',{1,2}}},...
                            {lb(1),lb(2)},{ub(1),ub(2)});
Constr2 = atlinconstraint(length(IMC1_md)+1,...
                          {{'alpha',{1}},...
                           {'alpha',{2}}},...
                          [0 0], [0 0], [1 1]);

[IMC1_mdb,penalty_mdb,dmin_mdb] = atmatch_mod(IMC1_md,Variab3,Constr2,...
                     TolAlphas,TolX,100,verboselevel,@fminsearch,twiss_uc);
if (strcmp(verbosef,'Y'))
       toc;
       fprintf('Alphas matched with penalty = %8.3e \n',sum(penalty_mdb.^2));
       disp(dmin_mdb);
end

if (strcmp(plotf,'Y'))
    lindata=atlinopt(IMC1_mdb,0.0,1:length(IMC1_mdb)+1, 'twiss_in', twiss_uc);
    PlotBetaDisp(lindata,'Dispersion+Alpha Match');
end
%
% if fit is not accepted, then try a joint fit of dispersion and alphas
%
IMC1_mdbj = IMC1_mdb;

% only checks for jont fits if nojointf is 'N'
if (not(strcmp(nojointf,'Y')))
    disjointfit = questdlg('Keep this fit ?');
else
    disjointfit = 'Yes';
end

if(not(strcmp(disjointfit,'Yes')))
%
% Alternative match of dispersion and alphas simultaneously
%
    if (strcmp(verbosef,'Y'))
        tic;
        fprintf('Joint matching of dispersion and alphas \n');
    end

    Variables   = [Variab1a,Variab1b,Variab2,Variab3];
    Constraints = [Constr1,Constr2];

    [IMC1_mdbj, penalty_mdbj, dmin_mdbj] = atmatch_mod(IMC1_T,Variables,...
                         Constraints,TolX, TolDisp,nit,...
                         verboselevel,@fminsearch,...
                         twiss_uc);
    if (strcmp(verbosef,'Y'))
        toc;
        fprintf('Alphas and dispersion matched with penalty = %8.3e \n', sum(penalty_mdbj.^2));
        disp(dmin_mdbj);
    end
    if (strcmp(plotf,'Y'))
        lindata=atlinopt(IMC1_mdbj,0.0,1:length(IMC1_mdbj)+1,... 
        'twiss_in', twiss_uc);
        PlotBetaDisp(lindata,'Dispersion+Alpha Joint Match');
    end
end

%twiss_0 = atlinopt(IMC1_mdbj, 0.0, length(IMC1_mdbj) + 1, 'twiss_in', twiss_uc);
%MC1 = flip(IMC1_mdbj);

%
% sets the new quad strengths to the whole achromat
%
DVsMC  = getDVs(4,IMC1_mdb,LatticeOptData);
DVsUC  = getDVs(3, UC_T, LatticeOptData);
ACHRO  = setDVs(2,ACHRO,LatticeOptData,DVsMC);
ACHRO  = setDVs(2,ACHRO,LatticeOptData,DVsUC);
DVs    = getDVs(2,ACHRO,LatticeOptData);

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

if (strcmp(plotf,'Y'))
    lindata_achr=atlinopt(ACHRO,0.0,1:length(ACHRO)+1);
    PlotBetaDisp(lindata_achr,'Full Achromat');
end

if (strcmp(verbosef,'Yes'))
    disp('Calculating optics parameters...');
    tic;
end
try
    pars=atsummary_fast(ACHRO,LatticeOptData.isdipole);
catch
    disp('Optics calculation not possible...');
    pars.circ=NaN;
    pars.naturalEmittance=NaN;
    pars.naturalEnergySpread=NaN;
    pars.tunes=[NaN NaN];
    pars.damping = [NaN NaN NaN];
    pars.chromaticity= [NaN NaN];
    pars.alphac = NaN;
    pars.Qx_ring = NaN;
    pars.Qy_ring = NaN;
    pars.beta0 = [NaN NaN];
    pars.etax0 = NaN;
end
if(strcmp(verbosef,'Y'))
    toc;
end
dCirc = (pars.circ-528.0)*1000;
dfRF  = dCirc/528.0*frf/1000;
unpkick = RbLength*Krb*yOffset;
RbShift = (dTheta*1000 - unpkick)/RbLength/Krb;

Ksixring = pars.chromaticity(1)*20;
Ksiyring = pars.chromaticity(2)*20;

if (strcmp(verbosef,'Y'))
    disp ('Calculating twiss...');
    tic;
end

if(strcmp(verbosef,'Y'))
    toc;
end
%
rp.inputs.Trb             = dTheta;
rp.inputs.disjointfit     = disjointfit;
rp.inputs.LatticeOptData  = LatticeOptData;
rp.inputs.uctunes         = uctunes;
rp.inputs.nit             = nit;
rp.inputs.TolTunes        = TolTunes;
rp.inputs.TolX            = TolX;
rp.inputs.TolDisp         = TolDisp;
rp.inputs.TolAlphas       = TolAlphas;
rp.inputs.lb              = lb;
rp.inputs.ub              = ub;

rp.outputs.ringpar          = pars;
rp.outputs.ringpar.dfRF     = dfRF;
rp.outputs.ringpar.yOffset  = yOffset;
rp.outputs.ringpar.unpkick  = unpkick;
rp.outputs.ringpar.RbShitf  = RbShift;
rp.outputs.ringpar.Ksixring = Ksixring;
rp.outputs.ringpar.Ksiyring = Ksiyring;

rp.outputs.DVs             = DVs;

if(strcmp(saveOPAf,'Y'))
   DVs=getDVs(2,ACHRO,LatticeOptData); 
   filein   = 'm4_Studies_InputfromMOGA_AT_Template_COMP.opa';
   fileout  = 'm4_Studies_InputfromMOGA_AT_COMP.opa';
   fileIDin = fopen(filein,'r');
   pathnameout = 'C:\Users\pedtav\Documents\Accelerators\NewMachineStudies\Lattices\OrbitShift\OPA\';
   fileIDout = fopen([pathnameout fileout],'w');
   fprintf(fileIDout, '{ Decision Variables from Matched Lattice }\n');
   for i=1:nvars
       fprintf(fileIDout, 'DV%1d = %8.5f ; \n', i, DVs(i));
   end
   fprintf(fileIDout, 'RBK_MOGA = %8.5f ;\n', dTheta);
    
   while 1
        tline = fgetl(fileIDin);
        if ~ischar(tline), break, end
        fprintf(fileIDout,'%s \n', tline);
    end
     fclose(fileIDin);
     fclose(fileIDout);
     fprintf('saved OPA file %s \n',fileout);
end       
           
end

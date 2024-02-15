function rbpars = RBParam_AT2(dTheta, nuc, varargin)
%RBParam 
%   Returns beam parameters for ring composed only of unit cells
%   with two reverse bends.
%   Input:  
%          dTheta = Reverse bend angle  [rad]
%          nuc = number of unit cells for calculation of RDTs
%          varargin{1} = targetchrom = [ksix, ksiy] target chromaticties
%          varargin{2} reverse bend quadrupolar strength - default is same
%          as in lattice without reverse bend [m**-2]
%
%   UCRB is a global cell array containing the lattice of one unit cell
%
%   Ouput: [dTheta*1000 Krb Circ dCirc dfRF 0.0 yOffset unpkick RbShift ...
%           emit Remit Jx alphac sigma_E Qx Qy Ksix Ksiy betax0 betay0 disp0 RDTs apx amx apy];
%   dTheta = Reverse bend angle [mrad]
%   Krb    = Reverse Bend Quadrupole strength [m**-2]
%   Circ   = Circumference[m]
%   dCirC  = change in circumference compare to zero RB kick case [mm]
%   dfRF   = change in RF frequency [Hz]
%   TotAngle = total bend angle [deg]
%   yOffset  = vertical orbvoit offset at the position of the Reverse Bends [mm]
%   unpkick  = kick from Reverse bend due to the vertical offset above [mrad]
%   RbShift  = additional horiuzonatl shift of Reverse bend in order to
%   achieve requierd total RB kick [mm]
%   emit     = emittance [pm rad]
%   Remit    = Ratio of emittance to unperturbed case [%]
%   Jx       = Horizontal damping partition
%   alphac   = momentum compaction [0.001]
%   sigma_E  = energy spread [0.001]
%   Qx, Qy   = Betatron Tunes
%   Ksix Ksiy = Chormaticities
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
%
global UCRB RUCRB;

frf=99.931E6;
RbLength = 0.15;

if (nargin>2)
    chroms = varargin{1};
    fitchrom='Yes';
else
    fitchrom='No';
end
    
if (nargin>3)
    Krb = varargin{2};
    max4_UpgradeStudies_20230305_AT2(nuc,'Yes',dTheta,Krb);
else
    max4_UpgradeStudies_20230305_AT2(nuc,'Yes',dTheta);
    Krb = 4.030076;
end
disp('Calculating orbit geometry...');
tic;
[x2d y2d a2d] = Survey2D(UCRB,1.5*pi/180);
yOffset = 6.509501134666746-y2d(round(length(y2d)/2))*1000;
toc;
if (strcmp(fitchrom,'Yes')==1)
    disp('Fitting Chromaticity...')
    tic
    try
        UCRB_C=atfitchrom(UCRB,chroms, 'sd', 'sfi'); % careful: a single iteration may not be enough
    catch
        disp('Chromaticity fitting not possible');
        UCRB_C=UCRB;
    end
    toc
else
    UCRB_C=UCRB;
end

% finds the sextupole strengths in the unit cell
%
I_sds  = find(atgetcells(UCRB_C,'FamName','sd'));
I_sfis = find(atgetcells(UCRB_C,'FamName','sfi'));
K_sds   = atgetfieldvalues(UCRB_C,I_sds(1),'PolynomB',{1,3});
K_sfis = atgetfieldvalues(UCRB_C,I_sfis(1),'PolynomB',{1,3});

% sets sextupole strengths for the ring composed of unit cells
%
I_sdsR = find(atgetcells(RUCRB,'FamName','sd'));
I_sfisR = find(atgetcells(RUCRB,'FamName','sfi'));
RUCRB_C=atsetfieldvalues(RUCRB,I_sdsR,'PolynomB',{1,3},K_sds);
RUCRB_C=atsetfieldvalues(RUCRB_C,I_sfisR,'PolynomB',{1,3},K_sfis);


disp('Calculating optics parameters...');
tic;
try
    pars=ringpara(UCRB_C);
catch
    disp('Optics calculation not possible...');
    pars.R=NaN;
    pars.emittx=NaN;
    pars.tunes=[NaN NaN];
    pars.dampingJ = [NaN NaN NaN];
    pars.chroms = [NaN NaN];
    pars.alphac = NaN;
    pars.sigma_E = NaN;
end
toc;
Circ = pars.R*2*pi*120;
dCirc = (Circ-336.000)*1000;
dfRF  = dCirc/336*frf/1000;
unpkick = RbLength*Krb*yOffset;
RbShift = (dTheta*1000 - unpkick)/RbLength/Krb;
emit = pars.emittx*1e12;
Remit = emit/336.9050205153787*100;
Qx = pars.tunes(1)*120;
Qy = pars.tunes(2)*120;
Jx  = pars.dampingJ(1);
Ksix = pars.chroms(1)*120;
Ksiy = pars.chroms(2)*120;
alphac = pars.alphac*1000;
sigma_E = pars.sigma_E*1000;
disp ('Calculating twiss');
tic;
try
    tw = twissring(UCRB_C,0,[1],'chrom');
    betax0 = tw.beta(1);
    betay0 = tw.beta(2);
    disp0  = tw.Dispersion(1);
catch
    disp('Twiss calculation not possible ...');
    betax0 = NaN;
    betay0 = NaN;
    disp0  = NaN;
end
toc;
%
%Resonance driving terms
%
 RDTs.h11001 = NaN;
 RDTs.h00111 = NaN;
 RDTs.h20001 = NaN;
 RDTs.h00201 = NaN;
 RDTs.h10002 = NaN;
 RDTs.h10010 = NaN;
 RDTs.h10100 = NaN;
 RDTs.h21000 = NaN;
 RDTs.h30000 = NaN;
 RDTs.h10110 = NaN;
 RDTs.h10020 = NaN;
 RDTs.h10200 = NaN;
 RDTs.h22000 = NaN;
 RDTs.h11110 = NaN;
 RDTs.h00220 = NaN;
 RDTs.h31000 = NaN;
 RDTs.h40000= NaN;
 RDTs.h20110 = NaN;
 RDTs.h11200 = NaN;
 RDTs.h20020 = NaN;
 RDTs.h20200 = NaN;
 RDTs.h00310 = NaN;
 RDTs.h00400 = NaN;
 RDTs.dnux_dJx = NaN;
 RDTs.dnux_dJy = NaN;
 RDTs.dnuy_dJy = NaN;
 
if (not(isnan(betax0)))
    disp('Calculating RDTs...');
    try
        RDTs   = computeRDT(RUCRB_C,1);
    catch
        disp('Not possible to calculate RDTs...');
    end
end

%
% Calculate dynamic aperture
% 
 apx=NaN;
 amx=NaN;
 apy=NaN;
if (not(isnan(betax0)))
    try
        [xx zz] = atdynap(UCRB_C, 100, 0.0, 0.02);
        xx=xx*1000;
        zz=zz*1000;
        apx=xx(1);
        amx=xx(length(xx));
        apy=zz((length(zz)+1)/2);
        figure;plot(xx,zz,'-o');
    catch
        disp('Error calculating Dynamic Aperture');
    end
end
    acpx = apx^2/betax0;
    acmx = amx^2/betax0;
    acy  = apy^2/betay0; 
    ac   = (acpx+acmx)*acy;
    
rbpars = [dTheta*1000 Krb Circ dCirc dfRF 0.0 yOffset unpkick RbShift ...
           emit Remit Jx alphac sigma_E Qx Qy Ksix Ksiy betax0 betay0 disp0 K_sds K_sfis ...
           abs(RDTs.h11001) ...
           abs(RDTs.h00111) ...
           abs(RDTs.h20001) ...
           abs(RDTs.h00201) ...
           abs(RDTs.h10002) ...
           abs(RDTs.h10010) ...
           abs(RDTs.h10100) ...
           abs(RDTs.h21000) ...
           abs(RDTs.h30000) ...
           abs(RDTs.h10110) ...
           abs(RDTs.h10020) ...
           abs(RDTs.h10200) ...
           abs(RDTs.h22000) ...
           abs(RDTs.h11110) ...
           abs(RDTs.h00220) ...
           abs(RDTs.h31000) ...
           abs(RDTs.h40000) ...
           abs(RDTs.h20110) ...
           abs(RDTs.h11200) ...
           abs(RDTs.h20020) ...
           abs(RDTs.h20200) ...
           abs(RDTs.h00310) ...
           abs(RDTs.h00400) ...
           RDTs.dnux_dJx ...
           RDTs.dnux_dJy ...
           RDTs.dnuy_dJy ...
           apx amx apy acpx acmx acy ac];
toc;
    

end

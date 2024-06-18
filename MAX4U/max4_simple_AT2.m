function varargout = max4_simple_AT2(varargin)
%MAX4_SIMPLE_AT2 Design date 20121107, branch 430, simplified version
%  Loads an Accelerator Toolbox lattice model of the MAX IV 3 GeV storage
%  ring into the cell array THERING, or the output.
%
%  Manually converted from TRACY-III. Source lattice file is
%  "m4-20121107-430-bare.lat". Field strengths in the dipole slices
%  correspond to the original theoretical ones.
%
%   The lattice includes:
%     - A sliced hard-edge approximation of the main dipole fringe-fields. 
%     - BPMs
%     - High bandwidth corrector magnets for FOFB (not present in source
%        lattice)
%     - Thin trim coils in the SDE, SFI, SFM and SFO sextupoles 
%     - Thin trim coils in the OYY octupoles
%     - Approximate models of the 3rd harmonic passive Landau cavities; see
%       comments below for limitations.
%     
%   The lattice does NOT include:
%     - Insertion devices
%     - Apertures at all restrictions (long straight chambers are included)
%     - Markers for diagnostic beamlines
%
%   Further comments:
%     - The lattice do not include or indicate which magnets are connected
%        in series. For this AT parameter groups are used.
%     - Support markers are NOT correctly placed and will be updated
%        at a later point in order to support misalignment/realignment
%        simulations.
%     - The harmonic cavities are approximated with a static longitudinal
%       kick. This is a valid approximation for the equilibrium particle
%       only. Note that the HC voltage after loading the lattice is zero if
%       default options are used.
%
%  In order to do realistic 6D tracking, classical radiation and the cavity
%  must be turned on using the cavityon, radiationon (AT1.3) commands,
%  alternatively atenable_6d (AT2.0).
%
%  Lattice designer Simon C. Leemann
%
%  Author: Magnus Sjöström, MAX IV laboratory (2014-11-11, rev. 2023-02-24).
%
%  See also cavityon, radiationon, fittune2, fitchrom2

%  Added scraper in achromat 11, David K. Olsson, 2017-06-09
%  Added second vertical scraper in achromat 1. Hopefully the positions
%  are correct, David K. Olsson, 2017-06-09
%  Added HC approximations, M. Sjöström, 2022-11-05


clear GLOBAL FAMLIST THERING GLOBVAL

str = dbstack;
fprintf('   Loading MAX IV 3 GeV ring lattice from file: %s\n',str(1).name);

E0 = 3000e6;
LatticeFile = str(1).name;

% RelGamma    = E0/(PhysConstant.electron_mass_energy_equivalent_in_MeV.value * 1e6);
% RelBeta     = sqrt(1 - RelGamma^-2);
RelBeta     = 1;                        % Not correct, but is implied to be an underlying assumption in AT2 tracking (see atparticle)
HarmNumber  = 176;                      % # of RF buckets in the ring
L0          = 528;
V0          = 1.8e6 / 6;                % Individual main cavity voltage, in unit V
HCvoltage   = 0;                        % Individual HC voltage, in unit V

% Vacuum chamber definition
chamberRadius = 0.011;      % Standard chamber inner radius
lsAxis = [0.018 0.004];     % Long straight Al chamber apertures, verified in the vacuum "picture book"
% bpmRadius = 0.0125;         % BPM chambers are steel, greater internal radius in order to be shadowed.

% BPM definitions
bpmhv     = atmonitor('BPM', 'IdentityPass','EApertures',[chamberRadius, chamberRadius]);
bpm_d     = atdrift('bpm_d', 0.025000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
hvpu      = cellcat(bpm_d, bpmhv, bpm_d);

% Extra BPM definition
% ehvpu     = atmarker('hvpue', 'IdentityPass');
% bpm_d     = atdrift('bpm_d', 0.025000,'DriftPass');
hvpue      = cellcat(bpm_d, bpmhv, bpm_d);


% Stripline definition
striphv     = atmarker('striphv', 'IdentityPass','EApertures',[chamberRadius, chamberRadius]);
strip_d     = atdrift('strip_d', 0.125000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
stripline   = cellcat(strip_d, striphv, strip_d);

% Scraper definition
scraperh    = ataperture('scraperh', [-1, 1, -1, 1],'AperturePass');
scraperv    = ataperture('scraperv', [-1, 1, -1, 1],'AperturePass');
scraper_d1  = atdrift('scraper_d', 0.118,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
scraper_d2  = atdrift('scraper_d', 0.100,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
scraper_d3  = atdrift('scraper_d', 0.118,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
% scraper_d1a = atdrift('str01064', 0.1064,'DriftPass');
% scraper_d2a = atdrift('str0035', 0.035,'DriftPass');
% scraper_d3a = atdrift('str00946', 0.0946,'DriftPass');
scraper     = cellcat(scraper_d1, scraperv, scraper_d2, scraperh, scraper_d3);
% scraper     = [scraper_d1a scraperh scraper_d2 scraperv scraper_d2a scraperv scraper_d3a};

% Add markers
% Girder/block start and end points use markers with FamName GS and GE.
% GirderStart = atmarker('GS','IdentityPass');
% GirderEnd   = atmarker('GE','IdentityPass');

GirderStart.M1 = atmarker('GS','IdentityPass','GirderType','M1');
GirderEnd.M1   = atmarker('GE','IdentityPass','GirderType','M1');
GirderStart.M2 = atmarker('GS','IdentityPass','GirderType','M2');
GirderEnd.M2   = atmarker('GE','IdentityPass','GirderType','M2');
GirderStart.U1 = atmarker('GS','IdentityPass','GirderType','U1');
GirderEnd.U1   = atmarker('GE','IdentityPass','GirderType','U1');
GirderStart.U2 = atmarker('GS','IdentityPass','GirderType','U2');
GirderEnd.U2   = atmarker('GE','IdentityPass','GirderType','U2');
GirderStart.U3 = atmarker('GS','IdentityPass','GirderType','U3');
GirderEnd.U3   = atmarker('GE','IdentityPass','GirderType','U3');
GirderStart.U4 = atmarker('GS','IdentityPass','GirderType','U4');
GirderEnd.U4   = atmarker('GE','IdentityPass','GirderType','U4');
GirderStart.U5 = atmarker('GS','IdentityPass','GirderType','U5');
GirderEnd.U5   = atmarker('GE','IdentityPass','GirderType','U5');


% Corrector definitions
corrh         =     atcorrector('corrh',0.0,[ 0 0 ],'CorrectorPass','iscorH','H');
corrv         =     atcorrector('corrv',0.0,[ 0 0 ],'CorrectorPass','iscorV','V');
corr_d        =     atdrift('corr_d', 0.025000,'DriftPass');
corr          =     cellcat(corr_d, corrv, corr_d, corr_d, corrh, corr_d);
corrs         =     cellcat(corr_d, corr_d, corr_d, corrh, corr_d);

corrh_ac    =   atcorrector('corrh_ac',0.0,[ 0 0 ],'CorrectorPass');
corrv_ac    =   atcorrector('corrv_ac',0.0,[ 0 0 ],'CorrectorPass');
corr_ac     =   cellcat(corrh_ac, corrv_ac);

% Absorber definitions

% RF definitions
CAV     = atrfcavity('CAV', 0.378, V0, HarmNumber*PhysConstant.speed_of_light_in_vacuum.value*RelBeta/L0, HarmNumber, E0, 'DriftPass','NumIntSteps',10,'EApertures',[chamberRadius, chamberRadius]);

% Landau passive cavity approximation
landau  = atbaselem('landau', 'StrMPoleSymplectic4Pass', ...
    'Length', 0.378, ...
    'K', 0, ...
    'PolynomA', [0,0,0,0], ...
    'PolynomB', [0,0,0,0], ...
    'T1', [0, 0, 0, 0, -HCvoltage / (2*E0), 0]', ...
    'T2', [0, 0, 0, 0, -HCvoltage / (2*E0), 0]', ...
    'NumIntSteps', 10, ...
    'EApertures',[chamberRadius, chamberRadius]);


% Long straight drift definitions
lstr0500  = atdrift('lstr0500', 0.500000,'DriftPass','EApertures',lsAxis);
lstr0321  = atdrift('lstr0321', 0.321000,'DriftPass','EApertures',lsAxis);

% Drift definitions
stra500  = atdrift('stra500', 0.500000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0500  = atdrift('str0500', 0.500000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0462  = atdrift('str0462', 0.462000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0420  = atdrift('str0420', 0.420000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
strx403  = atdrift('strx403', 0.403110,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0377  = atdrift('str0377', 0.377000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0362  = atdrift('str0362', 0.362000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0355  = atdrift('str0355', 0.355000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
% str0321  = atdrift('str0321', 0.321000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0302  = atdrift('str0302', 0.302000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0280  = atdrift('str0280', 0.280000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0270  = atdrift('str0270', 0.270000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0269  = atdrift('str0269', 0.269000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0252  = atdrift('str0252', 0.252000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
strx203  = atdrift('strx203', 0.203110,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0151  = atdrift('str0151', 0.151000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0150  = atdrift('str0150', 0.150000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0100  = atdrift('str0100', 0.100000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
strx093  = atdrift('strx093', 0.092670,'DriftPass','EApertures',[chamberRadius, chamberRadius]);    % NB! Length adjustment of 10 µm
str0075  = atdrift('str0075', 0.075000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
% str0050  = atdrift('str0050', 0.050000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0048  = atdrift('str0048', 0.048000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
stra025  = atdrift('stra025', 0.025000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0025  = atdrift('str0025', 0.025000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0021  = atdrift('str0021', 0.021000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0020  = atdrift('str0020', 0.020000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
strx013  = atdrift('strx013', 0.012500,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str0010  = atdrift('str0010', 0.010000,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
strx006  = atdrift('strx006', 0.006080,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str01836 = atdrift('str01836', 0.18360,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str15514 = atdrift('str15514', 1.55140,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str00594 = atdrift('str00594', 0.05940,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str03096 = atdrift('str03096', 0.30960,'DriftPass','EApertures',[chamberRadius, chamberRadius]);

% Straights added for achromat 11 scraper
str01925 = atdrift('str01925', 0.1925,'DriftPass','EApertures',[chamberRadius, chamberRadius]);
str03075 = atdrift('str03075', 0.3075,'DriftPass','EApertures',[chamberRadius, chamberRadius]);

% Septum definitions
septum  = cellcat(repmat({stra500},1,2), stra025);

% Injection kickers and pingers
ip      = atmarker('ip', 'IdentityPass','EApertures',[chamberRadius, chamberRadius]);

psmc    = atmultipole('psm'   , 0.0     , [0 0 0 0], [0 0 0*(0.9265e-3/6.057e-3/6.057e-3) 0], 'StrMPoleSymplectic4Pass','MaxOrder',3,'EApertures',[chamberRadius, chamberRadius]);
psm     = cellcat(str0150, psmc, str0150);

phc     = atcorrector('PMH3',0.0,[ 0 0 ],'IdentityPass');
pvc     = atcorrector('PMV3',0.0,[ 0 0 ],'IdentityPass');
ph      = cellcat(str0150, phc, str0150);
pv      = cellcat(str0075, pvc, str0075);

% Quadrupoles
% As some routines rely on the K value being set, equalize it to the
% PolynomB(2) value.
qf      = atmultipole('qf'    , 0.150000, [0 0 0 0], [0 4.03008 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[chamberRadius, chamberRadius]);
qfm   	= atmultipole('qfm'   , 0.150000, [0 0 0 0], [0 3.77399 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[chamberRadius, chamberRadius]);   
qfend 	= atmultipole('qfend' , 0.250000, [0 0 0 0], [0 3.65385 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[chamberRadius, chamberRadius]);  
qdend 	= atmultipole('qdend' , 0.250000, [0 0 0 0], [0 -2.50366 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[chamberRadius, chamberRadius]);   
qf.K      = qf.PolynomB(2);
qfm.K     = qfm.PolynomB(2);
qfend.K   = qfend.PolynomB(2);
qdend.K   = qdend.PolynomB(2);

% Sextupoles
sd_p    = atmultipole('sd'    , 2*0.050000, [0 0 0 0], [0 0 -116.62523 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[chamberRadius, chamberRadius]);
sdend_p	= atmultipole('sdend' , 2*0.050000, [0 0 0 0], [0 0 -170.000000 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[chamberRadius, chamberRadius]);   
sfm_p  	= atmultipole('sfm'   , 2*0.050000, [0 0 0 0], [0 0 170.000000 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[chamberRadius, chamberRadius]);   
sfo_p 	= atmultipole('sfo'   , 2*0.050000, [0 0 0 0], [0 0 174.000000 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[chamberRadius, chamberRadius]);  
sfi_p   = atmultipole('sfi'   , 2*0.050000, [0 0 0 0], [0 0 207.41203 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[chamberRadius, chamberRadius]); 
sd      = sd_p; %cellcat(sd_p, sdtrim, sd_p);
sdend   = sdend_p; %cellcat(sdend_p, sdetrim, sdend_p);
sfm     = sfm_p; %cellcat(sfm_p, sfmtrim, sfm_p);
sfo     = sfo_p; %cellcat(sfo_p, sfotrim, sfo_p);
sfi     = sfi_p; %cellcat(sfi_p, sfitrim, sfi_p);

% Octupoles - doublecheck whether strengths are integrated
oxxo  = atmultipole('oxxo' , 0.1, [0 0 0 0], [0 0 0 -1648.58], 'StrMPoleSymplectic4Pass','MaxOrder',3,'EApertures',[chamberRadius, chamberRadius]);  % Difference from OPA by an exact factor 10...
oxyo  = atmultipole('oxyo' , 0.1, [0 0 0 0], [0 0 0 3270.14], 'StrMPoleSymplectic4Pass','MaxOrder',3,'EApertures',[chamberRadius, chamberRadius]);   % Difference from OPA by an exact factor 10...
oyyo  = atmultipole('oyyo' , 0.1, [0 0 0 0], [0 0 0 -1420.22], 'StrMPoleSymplectic4Pass','MaxOrder',3,'EApertures',[chamberRadius, chamberRadius]);   % Difference from OPA by an exact factor 10...
% oxxo    = cellcat(str0050, oxxo, str0050);
% oxyo    = cellcat(str0050, oxyo, str0050);
% oyyo    = cellcat(str0050, oyyo, str0050);



% Bending magnet definitions
DIP.PolynomB = [ ...
                   0   -0.00013  0  0; ... DF5
                   0   0.01176  0  0; ... DF4
                   0  -0.55183  0  0; ... DF3
                   0  -0.86606  0  0; ... DF2
                   0  -0.86491  0  0; ... DF1
                   0  -0.86486  0  0; ... D0
                   0  -0.86486  0  0; ... D0 
                   0  -0.86491  0  0; ... DF1
                   0  -0.86606  0  0; ... DF2
                   0  -0.55183  0  0; ... DF3
                   0   0.01176  0  0; ... DF4
                   0   -0.00013  0  0]; %DF5
DIP.Length  = [ ...
    0.050000000000000; ...
    0.050000000000000; ...
    0.050000000000000; ...
    0.050000000000000; ...
    0.050000000000000; ...
    0.36190000000000; ...
    0.36190000000000; ...
    0.050000000000000; ...
    0.050000000000000; ...
    0.050000000000000; ...
    0.050000000000000; ...
    0.050000000000000];
DIP.BendingAngle = [ ...
    0.000089000000000; ... DF5
    0.001569000000000; ... DF4
    0.101861000000000; ... DF3
    0.151101000000000; ... DF2
    0.151199000000000; ... DF1
    1.094181000000000; ... D0
    1.094181000000000; ... D0
    0.151199000000000; ... DF1
    0.151101000000000; ... DF2
    0.101861000000000; ... DF3
    0.001569000000000; ... DF4
    0.000089000000000]; ... DF5

dipSlice = cell(1,size(DIP.PolynomB,1));               
for n = 1:size(DIP.PolynomB,1)
    dipSlice{n} = atsbend('dip', ...
        DIP.Length(n), ...
        DIP.BendingAngle(n)*pi/180, ...
        DIP.PolynomB(n,2), ...
        'BndMPoleSymplectic4Pass', ...
        'PolynomB', DIP.PolynomB(n,:), ...
        'MaxOrder',3, ...
        'EApertures',[chamberRadius, chamberRadius]);
end
   

DIPM.PolynomB = [ ...
                   0   0.00661   0  0; ...
                   0  -0.27143   0  0; ...
                   0  -0.42512   0  0; ...
                   0  -0.42605   0  0; ...
                   0  -0.58488   0  0; ...
                   0  -0.87035   0  0; ...
                   0  -0.87070   0  0; ...
                   0  -0.87075   0  0; ...
                   0  -0.87191   0  0; ...
                   0  -0.55556   0  0; ...
                   0   0.01184   0  0; ...
                   0  -0.00013   0  0];
DIPM.Length  = [ ...
   0.050000000000000; ...
   0.050000000000000; ...
   0.050000000000000; ...
   0.050000000000000; ...
   0.050000000000000; ...
   0.050000000000000; ...
   0.20420000000000; ...
   0.050000000000000; ...
   0.050000000000000; ...
   0.050000000000000; ...
   0.050000000000000; ...
   0.050000000000000];
DIPM.BendingAngle = [ ...   
   0.001070000000000; ...
   0.050729000000000; ...
   0.074672000000000; ...
   0.076248000000000; ...
   0.114983000000000; ...
   0.152049000000000; ...
   0.621695000000000; ...
   0.152220000000000; ...
   0.152122000000000; ...
   0.102549000000000; ...
   0.001579000000000; ...
   0.000090000000000];

dipmSlice = cell(1,size(DIPM.PolynomB,1));               
for n = 1:size(DIPM.PolynomB,1)
    dipmSlice{n} = atsbend('dipm', ...
        DIPM.Length(n), ...
        DIPM.BendingAngle(n)*pi/180, ...
        DIPM.PolynomB(n,2), ...
        'BndMPoleSymplectic4Pass', ...
        'PolynomB', DIPM.PolynomB(n,:), ...
        'MaxOrder',3, ...
        'EApertures',[chamberRadius, chamberRadius]);
end


% Gather dipole slices
dip  = dipSlice;
dipm = dipmSlice;

%% Define insertion devices

%  BALDER wiggler, built by SOLEIL. Total power 20.87 kW, Kmax=11.31,
%  Keff=8.99, N=38 plus 2 end pieces.
% Balder = atwiggler('BALDER', 2, 0.05, 0, E0, 'GWigSymplecticPass');

%% RING STRUCTURE 
%--------------------------------------------------------------------------


% Dipole magnet units
Dip     =   dip; %[dip, dip(end:-1:1)];
DipUC3  =   cellcat(dip{1:end/2}, dip{(end/2+1):end});

% Virtual units containing split quadrupoles with central sextupole
sqfm =      cellcat( qfm, str0075, sfm, strx013, hvpu, strx013, qfm );
sqfo =      cellcat( qf,  str0075, sfo, strx013, hvpu, strx013, qf );
sqfi =      cellcat( qf,  str0075, sfi, strx013, hvpu, strx013, qf );

% Construct straight sections
R3_01L_a        = cellcat( str01836, scraper, stripline, str15514 );
%R3_01L_b        = cellcat( str0500, str0377, septum, ip, str0098, str0321);
R3_01L_b        = cellcat( str0500, str0377, septum, ip, str00594, hvpue, str03096);
R3_02L_b        = cellcat( repmat({str0500},1,3), str0252, psm, str0269);
% DipKickerStr    = { repeat(str0500,4)', ph, str0021 );                        % Matches Simon's lattice, but not vacuum layout
R3_xxL          = cellcat( repmat({lstr0500},1,4), lstr0321);

ShortStr_RF     = cellcat( str0100, corr_ac, str0362, CAV, str0462 );           % Preliminary position for fast corrector, mounted on taper
ShortStr_HC     = cellcat( str0100, corr_ac, str0362, landau, str0462 );        % Preliminary position for fast corrector, mounted on taper

ShortStr_10_S2  = cellcat( str0355, str0302, pv, str0420, corr_ac, str0075 );   % Preliminary position for fast corrector, mounted on ion pump section
ShortStr_xx_S1  = cellcat( str0500, str0151, corr_ac, str0151, str0500 );       % Preliminary position for fast corrector, mounted on ion pump section
ShortStr_xx_S2  = cellcat( str0500, str0151, str0500, corr_ac, str0151 );       % Preliminary position for fast corrector, mounted on ion pump section

ShortStr_01_S1  = cellcat( str0270, ph, str0150, str0280, str0302, corr_ac );
ShortStr_02_S2  = cellcat( str0500, str0500, str0151, corr_ac, str0151 );       % Preliminary position for fast corrector, emittance monitor short straight
ShortStr_11_S2  = cellcat( str01925, scraperh, str03075, str0151, str0500, corr_ac, str0151);
ShortStr_20_S1  = cellcat( str0500, str0500, str0151, str0151 );                % Missing AC corrector for fast corrector, emittance monitor short straight

% Construct magnet sections
uc1     = cellcat( GirderStart.U1, sqfm, str0100, corr, strx203, sd, str0010, Dip, str0010, sd, GirderEnd.U1, strx403 );
uc2     = cellcat( GirderStart.U2, sqfo, str0100, corr, strx203, sd, str0010, Dip, str0010, sd, GirderEnd.U2, strx403 );
uc3     = cellcat( GirderStart.U3, sqfi, str0100, corr, strx203, sd, str0010, DipUC3, str0010, sd, ...
    strx203, corr, str0100, sqfi(end:-1:1), GirderEnd.U3 );
uc4     = cellcat( strx403, GirderStart.U4, sd, str0010, Dip(end:-1:1), str0010, sd, strx203, corr, str0100, sqfo(end:-1:1), GirderEnd.U4 );
uc5     = cellcat( strx403, GirderStart.U5, sd, str0010, Dip(end:-1:1), str0010, sd, strx203, corr, str0100, sqfm(end:-1:1), GirderEnd.U5 );


mc1a1   = cellcat( GirderStart.M1, hvpu, str0010, str0048, corr, str0021, oxxo, str0025, qfend, str0025, oxyo, str0100, ...
            qdend, strx006, dipm, oyyo, strx093, corrs, hvpu, str0020, sdend, GirderEnd.M1 );           % Missing AC corrector due to injection channel
mc2a1   = cellcat( GirderEnd.M2, hvpu, str0010, str0048, corr, str0021, oxxo, str0025, qfend, str0025, oxyo, str0100, ...
            qdend, strx006, dipm, oyyo, strx093, corr, hvpu, str0020, sdend, GirderStart.M2 );            % Missing AC corrector due to injection channel


mc1     = cellcat( GirderStart.M1, hvpu, str0010, corr_ac, str0048, corr, str0021, oxxo, str0025, qfend, str0025, oxyo, str0100, ...
            qdend, strx006, dipm, oyyo, strx093, corrs, hvpu, str0020, sdend, GirderEnd.M1 );
mc2     = cellcat( GirderEnd.M2, hvpu, str0010, corr_ac, str0048, corr, str0021, oxxo, str0025, qfend, str0025, oxyo, str0100, ...
            qdend, strx006, dipm, oyyo, strx093, corr, hvpu, str0020, sdend, GirderStart.M2 );

        
achr    = cellcat( R3_xxL, mc1, ShortStr_xx_S1, uc1, uc2, uc3, uc4, uc5, ShortStr_xx_S2, mc2(end:-1:1), R3_xxL(end:-1:1) );
achr1   = cellcat( R3_01L_b, mc1a1, ShortStr_01_S1, uc1, uc2, uc3, uc4, uc5, ShortStr_RF, mc2a1(end:-1:1), R3_xxL(end:-1:1) );
achr2   = cellcat( R3_02L_b, mc1, ShortStr_xx_S1, uc1, uc2, uc3, uc4, uc5, ShortStr_02_S2, mc2(end:-1:1), R3_xxL(end:-1:1) );
achr10  = cellcat( R3_xxL, mc1, ShortStr_xx_S1, uc1, uc2, uc3, uc4, uc5, ShortStr_10_S2, mc2(end:-1:1), R3_xxL(end:-1:1) );
achr11  = cellcat( R3_xxL, mc1, ShortStr_xx_S1, uc1, uc2, uc3, uc4, uc5, ShortStr_11_S2, mc2(end:-1:1), R3_xxL(end:-1:1) );
achr13  = cellcat( R3_xxL, mc1, ShortStr_xx_S1, uc1, uc2, uc3, uc4, uc5, ShortStr_HC, mc2(end:-1:1), R3_xxL(end:-1:1) );
achr14  = cellcat( R3_xxL, mc1, ShortStr_xx_S1, uc1, uc2, uc3, uc4, uc5, ShortStr_HC, mc2(end:-1:1), R3_xxL(end:-1:1) );
achr15  = cellcat( R3_xxL, mc1, ShortStr_xx_S1, uc1, uc2, uc3, uc4, uc5, ShortStr_HC, mc2(end:-1:1), R3_xxL(end:-1:1) );
achr16  = cellcat( R3_xxL, mc1, ShortStr_xx_S1, uc1, uc2, uc3, uc4, uc5, ShortStr_RF, mc2(end:-1:1), R3_xxL(end:-1:1) );
achr17  = cellcat( R3_xxL, mc1, ShortStr_xx_S1, uc1, uc2, uc3, uc4, uc5, ShortStr_RF, mc2(end:-1:1), R3_xxL(end:-1:1) );
achr18  = cellcat( R3_xxL, mc1, ShortStr_xx_S1, uc1, uc2, uc3, uc4, uc5, ShortStr_RF, mc2(end:-1:1), R3_xxL(end:-1:1) );
achr19  = cellcat( R3_xxL, mc1, ShortStr_xx_S1, uc1, uc2, uc3, uc4, uc5, ShortStr_RF, mc2(end:-1:1), R3_xxL(end:-1:1) );
achr20  = cellcat( R3_xxL, mc1, ShortStr_20_S1, uc1, uc2, uc3, uc4, uc5, ShortStr_RF, mc2(end:-1:1), R3_01L_a );

ring    = cellcat( achr1, atmarker('AchrEnd'), ...
            achr2, ...
            repmat(achr,1,7), ...
            achr10, ...
            achr11, ...
            achr, ...
            achr13, ...
            achr14, ...
            achr15, ...
            achr16, ...
            achr17, ...
            achr18, ...
            achr19, ...
            achr20 );

% For AT versions < 2.0 set Energy, MaxOrder and NumIntSteps globally.
% While AT 2.0 can use this information it is not sufficient to avoid
% slowdowns. In addition, set the customary global variables
for i=1:length(ring)
    ring{i}.Energy = E0;
    ring{i}.NumIntSteps = 20; 
    ring{i}.MaxOrder = 3;
end
global GLOBVAL
GLOBVAL.E0 = E0;
GLOBVAL.LatticeFile = LatticeFile;

% In order to avoid issues with orbit correction routines, which all seem
% to assume that we have as many vertical correctors as we do BPMs, we need
% to add weight information to the BPMs (and make sure correction routines
% make use of it).
I = findcells(ring,'FamName','BPM');
ring{I(1)}.Weight = [0 0]; I = I(2:end);
for n = 1:numel(I)
    if mod(n-2,10) == 0
        ring{I(n)}.Weight = [1 1e-3];
    else
        ring{I(n)}.Weight = [1 1];
    end
end

% Transpose the ring variable, as that seems to be what AT2 expects. Not
% that this won't risk causing problems with backward compatibility...
ring = ring(:); 

% If AT version is above 2.0, set the global properties. Required to avoid
% almost an order of magnitude slowdown due to look-ups and searches done
% routinely before tracking calculations
ring = atSetRingProperties(ring, ...
    'name',LatticeFile, ...
    'Energy',E0, ...                                % AT2 energy unit is eV
    'HarmNumber',HarmNumber, ...
    'rf_voltage', V0*numel(findcells(ring,'FamName','cav')), ...
    'particle','relativistic', ...                  % AT tracking assumes beta = 1, which this will result in
    'rf_frequency','nominal', ...                   % This recalculates the frequency based on the new rel. beta
    'cavpts',findcells(ring,'FamName','cav')+1);    % +1 required as the information is added as the first element in the cell array


% Calculate and print selected lattice properties
k = findspos(ring,length(ring)+1);
[RD, ~] = atlinopt(ring,0,1:numel(ring)+1);
nu = RD(end).mu/(2*pi);
[~, ksi] = tunechrom(ring,'get_chrom');
disp(['    R3 orbit length:   ' num2str(k,'%11.6g') ' m']);
disp(['    R3 betatron tunes: ' num2str(nu(1),'%3.6g') ' / ' num2str(nu(2),'%3.6g') '  (H/V)']);
disp(['    R3 chromaticities: ' num2str(ksi(1),'%3.6g') ' / ' num2str(ksi(2),'%3.6g') '  (H/V)']);

% If an output is expected, do not update THERING and instead send the
% lattice to the output
if nargout > 0
    varargout{1} = ring;
else
    global THERING
    THERING = ring;
    % Set up global parameters in the calling workspace, then exit
    evalin('caller','global THERING GLOBVAL ParamGroupDIP ParamGroupDIPM');
end
disp('    Loading completed...');



%% Helper functions

function catArray = cellcat(varargin)

N = 0;
for i = 1:numel(varargin)
    N = numel(varargin{i}) + N;
end
catArray = cell(1,N);

N = 1;
for i = 1:numel(varargin)
    if ~iscell(varargin{i})
        catArray{N} = varargin{i};
        N = N+1;
    else
        for j = 1:numel(varargin{i})
            catArray{N} = varargin{i}{j};
            N = N+1;
        end
    end
end
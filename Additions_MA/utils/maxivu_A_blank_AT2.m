function varargout = medmax_7BA_2_1_1_AT2_simple(varargin)
%MAX4_20121107_430E Design date 20121107, branch 430. Additional components added
%  Loads an Accelerator Toolbox lattice model of the MAX IV 3 GeV storage
%  ring into the cell array THERING.
%
%  Manually converted from TRACY-III. Source lattice file is
%  "m4-20121107-430-bare.lat". Field strengths in the dipole slices have
%  been changed to match the field profile determined in factory Hall probe
%  maps, after which the tunes and chromaticities were corrected.
%  Various elements such as scrapers, thick cavities, et.c. have been added
%  at their final locations.
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

global FAMLIST THERING GLOBVAL 
clear GLOBAL FAMLIST THERING GLOBVAL

PhysConst = PC;        %Load physical constants

str = dbstack;
fprintf('   Loading MAX IV 3 GeV ring lattice from file: %s\n',str(1).name);

E0 = 3000e6;
LatticeFile = str(1).name;

RelGamma    = E0/PhysConst.me_eV;
RelBeta     = sqrt(1 - RelGamma^-2);
% Ekin        = E0 - PhysConst.me_eV;
HarmNumber  = 176;                      % # of RF buckets in the ring
L0          = 528;
V0          = 1.8e6 / 6;                % Individual main cavity voltage, in unit V
HCvoltage   = 0;                        % Individual HC voltage, in unit V

% Vacuum chamber definition
CHAMBERs=   ataperture('CHAMBER', [-1, 1, -1, 1], 'EAperturePass');
CHAMBER.Axes = [0.035 0.0055];                     %The given values are "half-gaps"


% Add markers
% ....

% Corrector definitions
corrh         =     atcorrector('corrh',0.0,[ 0 0 ],'CorrectorPass');
corrv         =     atcorrector('corrv',0.0,[ 0 0 ],'CorrectorPass');
corr_d        =     atdrift('corr_d', 0.025000,'DriftPass');
corr          =     cellcat(corr_d, corrv, corr_d, corr_d, corrh, corr_d);
corrs         =     cellcat(corr_d, corr_d, corr_d, corrh, corr_d);

corrh_ac    =   atcorrector('corrh_ac',0.0,[ 0 0 ],'CorrectorPass');
corrv_ac    =   atcorrector('corrv_ac',0.0,[ 0 0 ],'CorrectorPass');
corr_ac     =   cellcat(corrh_ac, corrv_ac);

% Absorber definitions

% RF definitions
CAV     = atrfcavity('CAV', 0.378, V0, HarmNumber*PhysConst.c*RelBeta/L0, HarmNumber, E0, 'DriftPass','NumIntSteps',10);

% Landau passive cavity approximation
landau  = atbaselem('landau', 'QuadLinearPass', ...
    'Length', 0.378, ...
    'K', 0, ...
    'PolynomA', [0,0,0,0], ...
    'PolynomB', [0,0,0,0], ...
    'T1', [0, 0, 0, 0, -HCvoltage / (2*E0), 0], ...
    'T2', [0, 0, 0, 0, -HCvoltage / (2*E0), 0], ...
    'NumIntSteps', 10);


% Straight Section definitions
stra500  = atdrift('stra500', 0.500000,'DriftPass');
str0500  = atdrift('str0500', 0.500000,'DriftPass');
str0462  = atdrift('str0462', 0.462000,'DriftPass');
str0420  = atdrift('str0420', 0.420000,'DriftPass');
strx403  = atdrift('strx403', 0.403110,'DriftPass');
str0377  = atdrift('str0377', 0.377000,'DriftPass');
str0362  = atdrift('str0362', 0.362000,'DriftPass');
str0355  = atdrift('str0355', 0.355000,'DriftPass');
str0321  = atdrift('str0321', 0.321000,'DriftPass');
str0302  = atdrift('str0302', 0.302000,'DriftPass');

str0280  = atdrift('str0280', 0.280000,'DriftPass');
str0270  = atdrift('str0270', 0.270000,'DriftPass');
str0269  = atdrift('str0269', 0.269000,'DriftPass');
str0252  = atdrift('str0252', 0.252000,'DriftPass');
strx203  = atdrift('strx203', 0.203110,'DriftPass');
str0151  = atdrift('str0151', 0.151000,'DriftPass');
str0150  = atdrift('str0150', 0.150000,'DriftPass');
str0100  = atdrift('str0100', 0.100000,'DriftPass');
strx093  = atdrift('strx093', 0.092680,'DriftPass');
str0098  = atdrift('str0098', 0.098000,'DriftPass');
str0075  = atdrift('str0075', 0.075000,'DriftPass');
str0058  = atdrift('str0058', 0.058000,'DriftPass');
str0048  = atdrift('str0048', 0.048000,'DriftPass');
stra025  = atdrift('stra025', 0.025000,'DriftPass');
str0050  = atdrift('str0050', 0.050000,'DriftPass');
str0025  = atdrift('str0025', 0.025000,'DriftPass');
str0021  = atdrift('str0021', 0.021000,'DriftPass');
str0020  = atdrift('str0020', 0.020000,'DriftPass');
strx013  = atdrift('strx013', 0.012500,'DriftPass');
str0010  = atdrift('str0010', 0.010000,'DriftPass');
strx006  = atdrift('strx006', 0.006080,'DriftPass');
str01836 = atdrift('str01836', 0.18360,'DriftPass');
str15514 = atdrift('str15514', 1.55140,'DriftPass');
str00594 = atdrift('str00594', 0.05940,'DriftPass');
str03096 = atdrift('str03096', 0.30960,'DriftPass');

% Straights added for achromat 11 scraper
str01925 = atdrift('str01925', 0.1925,'DriftPass');
str03075 = atdrift('str03075', 0.3075,'DriftPass');


% Quadrupoles
% As some routines rely on the K value being set, equalize it to the
% PolynomB(2) value.  QF = 4.030928383463235
QF      = atmultipole('QF'    , 0.150000, [0 0 0 0], [0  4.030076 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2);
QFI     = atmultipole('QFI'   , 0.150000, [0 0 0 0], [0  4.030076 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2);
QFO     = atmultipole('QFO'   , 0.150000, [0 0 0 0], [0  4.030076 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2);
QFM   	= atmultipole('QFM'   , 0.150000, [0 0 0 0], [0  3.773995 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2);   
QFEND 	= atmultipole('QFEND' , 0.250000, [0 0 0 0], [0  3.653849 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2);  
QDEND 	= atmultipole('QDEND' , 0.250000, [0 0 0 0], [0 -2.503663 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2);   
QF.K      = QF.PolynomB(2);
QFI.K     = QF.PolynomB(2);
QFO.K     = QF.PolynomB(2);
QFM.K     = QFM.PolynomB(2);
QFEND.K   = QFEND.PolynomB(2);
QDEND.K   = QDEND.PolynomB(2);

%
% SEXTs
%
SD    = atmultipole('SD'    , 0.10000, [0 0 0 0], [0 0 -116.625229 0], 'StrMPoleSymplectic4Pass','MaxOrder',2);
SDEND = atmultipole('SDEND' , 0.10000, [0 0 0 0], [0 0 -170.000000 0], 'StrMPoleSymplectic4Pass','MaxOrder',2);   
SFM   = atmultipole('SFM'   , 0.10000, [0 0 0 0], [0 0  170.000000 0], 'StrMPoleSymplectic4Pass','MaxOrder',2);   
SFO   = atmultipole('SFO'   , 0.10000, [0 0 0 0], [0 0  174.000000 0], 'StrMPoleSymplectic4Pass','MaxOrder',2);  
SFI   = atmultipole('SFI'   , 0.10000, [0 0 0 0], [0 0  207.412038 0], 'StrMPoleSymplectic4Pass','MaxOrder',2); 
          
%
% OCTs
%
OXXO  = atmultipole('OXXO' , 0.10000, [0 0 0 0], [0 0 0 -1648.58], 'StrMPoleSymplectic4Pass','MaxOrder',3);  
OXYO  = atmultipole('OXYO' , 0.10000, [0 0 0 0], [0 0 0 3270.14],  'StrMPoleSymplectic4Pass','MaxOrder',3); 
OYYO  = atmultipole('OYYO' , 0.10000, [0 0 0 0], [0 0 0 -1420.22], 'StrMPoleSymplectic4Pass','MaxOrder',3); 

%
% BENDS
% 

D0  = atsbend('DIP', 0.361890, 1.094181*pi/180, -0.864858, 'BndMPoleSymplectic4Pass');
DF1 = atsbend('DIP', 0.050000, 0.151199*pi/180, -0.864908, 'BndMPoleSymplectic4Pass');
DF2 = atsbend('DIP', 0.050000, 0.151101*pi/180, -0.866059, 'BndMPoleSymplectic4Pass');
DF3 = atsbend('DIP', 0.050000, 0.101861*pi/180, -0.551829, 'BndMPoleSymplectic4Pass');
DF4 = atsbend('DIP', 0.050000, 0.001569*pi/180,  0.011759, 'BndMPoleSymplectic4Pass');
DF5 = atsbend('DIP', 0.050000, 0.000089*pi/180, -0.000128, 'BndMPoleSymplectic4Pass');

DS6 = atsbend('DIPm', 0.050000, 0.001070*pi/180,  0.006608, 'BndMPoleSymplectic4Pass');
DS5 = atsbend('DIPm', 0.050000, 0.050729*pi/180, -0.271428, 'BndMPoleSymplectic4Pass');
DS4 = atsbend('DIPm', 0.050000, 0.074672*pi/180, -0.425119, 'BndMPoleSymplectic4Pass');
DS3 = atsbend('DIPm', 0.050000, 0.076248*pi/180, -0.426048, 'BndMPoleSymplectic4Pass');
DS2 = atsbend('DIPm', 0.050000, 0.114983*pi/180, -0.584884, 'BndMPoleSymplectic4Pass');
DS1 = atsbend('DIPm', 0.050000, 0.152049*pi/180, -0.870351, 'BndMPoleSymplectic4Pass');
DS0 = atsbend('DIPm', 0.204240, 0.621695*pi/180, -0.870701, 'BndMPoleSymplectic4Pass');

          
DM1 = atsbend('DIPm', 0.050000, 0.152220*pi/180, -0.870751, 'BndMPoleSymplectic4Pass');
DM2 = atsbend('DIPm', 0.050000, 0.152122*pi/180, -0.871910, 'BndMPoleSymplectic4Pass');
DM3 = atsbend('DIPm', 0.050000, 0.102549*pi/180, -0.555557, 'BndMPoleSymplectic4Pass');
DM4 = atsbend('DIPm', 0.050000, 0.001579*pi/180,  0.011839, 'BndMPoleSymplectic4Pass');
DM5 = atsbend('DIPm', 0.050000, 0.000090*pi/180, -0.000129, 'BndMPoleSymplectic4Pass');

DIP  = cellcat(DF5, DF4, DF3, DF2, DF1, D0);
DIPm = cellcat(DS6, DS5, DS4, DS3, DS2, DS1, DS0, DM1, DM2, DM3, DM4, DM5);

ShortStr_RF     = cellcat( str0100, corr_ac, str0362, CAV, str0462 );           % Preliminary position for fast corrector, mounted on taper


BPM     = cellcat( str0025, str0025);
CORR    = cellcat( str0025, str0050, str0025);
CORRS   = cellcat( str0025, str0075);
PSM     = cellcat( str0150, str0150);
KI      = cellcat( str0150, str0150);
PV      = cellcat( str0075, str0075);
SEPT    = cellcat( stra500, stra500, stra025);

SQFM = cellcat(QFM, str0075, SFM, strx013, BPM, strx013, QFM, str0100, CORR);
SQFO = cellcat(QFO, str0075, SFO, strx013, BPM, strx013, QFO, str0100, CORR);
SQFI = cellcat(QFI, str0075, SFI, strx013, BPM, strx013, QFI, str0100, CORR);

DIPUC = cellcat(SD, str0010, DIP, flip(DIP), str0010, SD);

% OXX   = cellcat(str0050, OXXO, str0050);
% OXY   = cellcat(str0050, OXYO, str0050);
% OYY   = cellcat(str0050, OYYO, str0050);



LS      = cellcat( repmat(cellcat(str0500),1,4), str0321);
LS1A    = cellcat( str0500, str0377, SEPT, str0098, str0321);
LS2A    = cellcat( repmat(cellcat(str0500),1,3), str0252, PSM, str0269);
SS      = cellcat( repmat(cellcat(str0500),1,2), str0302);
SS1A    = cellcat( str0270, KI, str0150, str0280, str0302);
SS10B   = cellcat( str0420, str0075, PV, str0355, str0302);
UC1     = cellcat( SQFM, strx203, DIPUC, strx403);
UC2     = cellcat( SQFO, strx203, DIPUC, strx403);
UC3     = cellcat( SQFI, strx203, DIPUC, strx203, flip(SQFI));
UC4     = cellcat( strx403, DIPUC, strx203, flip(SQFO));
UC5     = cellcat( strx403, DIPUC, strx203, flip(SQFM));
MC1     = cellcat( BPM, str0058, flip(CORR), str0021, OXXO, str0025, QFEND, str0025, OXYO,...
          str0100, QDEND, strx006, DIPm, OYYO, strx093, CORRS, BPM, str0020,...
          SDEND);
MC2     = cellcat( BPM, str0058, flip(CORR), str0021, OXXO, str0025, QFEND, str0025, OXYO,...
          str0100, QDEND, strx006, DIPm, OYYO, strx093, CORR, BPM, str0020,...
          SDEND);
% ACHR1   = cellcat( LS1A, MC1, SS1A, UC1, UC2, UC3, UC4, UC5, flip(SS), flip(MC2), flip(LS));
% ACHR2   = cellcat( LS2A, MC1, SS, UC1, UC2, UC3, UC4, UC5, flip(SS), flip(MC2), flip(LS)) ;
% ACHR10  = cellcat( LS, MC1, SS, UC1, UC2, UC3, UC4, UC5, flip(SS10B), flip(MC2), flip(LS));
achrRF  = cellcat( LS, MC1, SS, UC1, UC2, UC3, UC4, UC5, ShortStr_RF, flip(MC2), flip(LS) );


achr    = cellcat( LS, MC1, SS, UC1, UC2, UC3, UC4, UC5, flip(SS), flip(MC2), flip(LS));
ring    = cellcat( achrRF, repmat(achr,1,15), repmat(achr,1,4)); % , repmat(achr,1, 5));
%RINGINJ = cellcat( ACHR1, ACHR2, 7*ACHR, ACHR10, 10*ACHR);

% RING STRUCTURE 
%--------------------------------------------------------------------------



% For AT versions < 2.0 set Energy, MaxOrder and NumIntSteps globally.
% While AT 2.0 can use this information it is not sufficient to avoid
% slowdowns. In addition, set the customary global variables
for i=1:length(ring)
    ring{i}.Energy = E0;
    ring{i}.NumIntSteps = 10; 
    ring{i}.MaxOrder = 3;
end
global GLOBVAL
GLOBVAL.E0 = E0;
GLOBVAL.LatticeFile = LatticeFile;


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
    'rf_frequency','nominal', ...
    'rf_voltage', V0*numel(findcells(ring,'FamName','cav')), ...
    'particle','relativistic', ...                  % AT tracking assumes beta = 1
    'cavpts',findcells(ring,'FamName','cav')+1);    % +1 required as the information is added as the first element in the cell array

% If an output is expected, do not update THERING and instead send the
% lattice to the output
if nargout > 0
    varargout{1} = ring;
else
    global THERING
    THERING = ring;
end

% Set up global parameters in the calling workspace, then exit
k = findspos(ring,length(ring)+1);
[RD, ~] = atlinopt(ring,0,1:numel(ring)+1);
nu = RD(end).mu/(2*pi);
[~, ksi] = tunechrom(ring,'get_chrom');
disp(['    R3 orbit length:   ' num2str(k,'%11.6g') ' m']);
disp(['    R3 betatron tunes: ' num2str(nu(1),'%3.6g') ' / ' num2str(nu(2),'%3.6g') '  (H/V)']);
disp(['    R3 chromaticities: ' num2str(ksi(1),'%3.6g') ' / ' num2str(ksi(2),'%3.6g') '  (H/V)']);
evalin('caller','global THERING GLOBVAL ParamGroupDIP ParamGroupDIPM');
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
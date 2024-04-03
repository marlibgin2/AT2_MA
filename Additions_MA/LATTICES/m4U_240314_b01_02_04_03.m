function varargout = m4U_240314_b01_02_04_03(varargin)

%
% lattice corresponding to the best so far 'B' case with
% with RB in U3
% stems from m4U_b1_2_1 where the number of power supplies
% has been incremented to control more gradients
%
%  
%

%%%global FAMLIST THERING GLOBVAL  
%%%clear GLOBAL FAMLIST THERING GLOBVAL

% obsolete, removed 15/2/2024 PhysConst = PC;        %Load physical constants
me_eV    = PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e6;
c_light  = PhysConstant.speed_of_light_in_vacuum.value;

str = dbstack;
fprintf('   Loading MAX IV 3 GeV ring lattice from file: %s\n',str(1).name);

E0 = 3000e6;
LatticeFile = str(1).name;

RelGamma    = E0/me_eV;
RelBeta     = sqrt(1 - RelGamma^-2);
% Ekin        = E0 - PhysConst.me_eV;
HarmNumber  = 176;                      % # of RF buckets in the ring
L0          = 528;
V0          = 1.8e6; % / 6;                % Individual main cavity voltage, in unit V
HCvoltage   = 0;                        % Individual HC voltage, in unit V

% Vacuum chamber definition
CHAMBERs     =   ataperture('CHAMBER', [-1, 1, -1, 1], 'EAperturePass');
CHAMBER.Axes =   [0.035 0.0055];                     %The given values are "half-gaps"


% Add markers
% ....

% Girder definition (start-end)
GS = atmarker('GS', 'IdentityPass'); % girder start marker
GE = atmarker('GE', 'IdentityPass'); % girder end marker 
GRD  = atmarker('GRD',  'IdentityPass'); % girder marker (either start or end)

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
%CAV     = atrfcavity('CAV', 0.378, V0, HarmNumber*PhysConst.c*RelBeta/L0, HarmNumber, E0, 'DriftPass','NumIntSteps',10);
CAV     = atrfcavity('CAV', 0.378, V0, HarmNumber*c_light*RelBeta/L0, HarmNumber, E0, 'DriftPass','NumIntSteps',10);

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
Rx = 0.011; Ry = 0.011; 
stra500  = atdrift('stra500', 0.500000,'DriftPass','EApertures',[Rx, Ry]);
str0500  = atdrift('str0500', 0.500000,'DriftPass','EApertures',[Rx, Ry]);
str0462  = atdrift('str0462', 0.462000,'DriftPass','EApertures',[Rx, Ry]);
str0420  = atdrift('str0420', 0.420000,'DriftPass','EApertures',[Rx, Ry]);
strx403  = atdrift('strx403', 0.403110,'DriftPass','EApertures',[Rx, Ry]);
str0377  = atdrift('str0377', 0.377000,'DriftPass','EApertures',[Rx, Ry]);
str0362  = atdrift('str0362', 0.362000,'DriftPass','EApertures',[Rx, Ry]);
str0355  = atdrift('str0355', 0.355000,'DriftPass','EApertures',[Rx, Ry]);
str0321  = atdrift('str0321', 0.321000,'DriftPass','EApertures',[Rx, Ry]);
str0302  = atdrift('str0302', 0.302000,'DriftPass','EApertures',[Rx, Ry]);

str0280  = atdrift('str0280', 0.280000,'DriftPass','EApertures',[Rx, Ry]);
str0270  = atdrift('str0270', 0.270000,'DriftPass','EApertures',[Rx, Ry]);
str0269  = atdrift('str0269', 0.269000,'DriftPass','EApertures',[Rx, Ry]);
str0252  = atdrift('str0252', 0.252000,'DriftPass','EApertures',[Rx, Ry]);
strx203  = atdrift('strx203', 0.203110,'DriftPass','EApertures',[Rx, Ry]);
str0151  = atdrift('str0151', 0.151000,'DriftPass','EApertures',[Rx, Ry]);
str0150  = atdrift('str0150', 0.150000,'DriftPass','EApertures',[Rx, Ry]);
str0100  = atdrift('str0100', 0.100000,'DriftPass','EApertures',[Rx, Ry]);
strx093  = atdrift('strx093', 0.092680,'DriftPass','EApertures',[Rx, Ry]);
str0098  = atdrift('str0098', 0.098000,'DriftPass','EApertures',[Rx, Ry]);
str0075  = atdrift('str0075', 0.075000,'DriftPass','EApertures',[Rx, Ry]);
str0058  = atdrift('str0058', 0.058000,'DriftPass','EApertures',[Rx, Ry]);
str0048  = atdrift('str0048', 0.048000,'DriftPass','EApertures',[Rx, Ry]);
stra025  = atdrift('stra025', 0.025000,'DriftPass','EApertures',[Rx, Ry]);
str0050  = atdrift('str0050', 0.050000,'DriftPass','EApertures',[Rx, Ry]);
str0025  = atdrift('str0025', 0.025000,'DriftPass','EApertures',[Rx, Ry]);
str0021  = atdrift('str0021', 0.021000,'DriftPass','EApertures',[Rx, Ry]);
str0020  = atdrift('str0020', 0.020000,'DriftPass','EApertures',[Rx, Ry]);
strx013  = atdrift('strx013', 0.012500,'DriftPass','EApertures',[Rx, Ry]);
str0010  = atdrift('str0010', 0.010000,'DriftPass','EApertures',[Rx, Ry]);
strx006  = atdrift('strx006', 0.006080,'DriftPass','EApertures',[Rx, Ry]);
str01836 = atdrift('str01836', 0.18360,'DriftPass','EApertures',[Rx, Ry]);
str15514 = atdrift('str15514', 1.55140,'DriftPass','EApertures',[Rx, Ry]);
str00594 = atdrift('str00594', 0.05940,'DriftPass','EApertures',[Rx, Ry]);
str03096 = atdrift('str03096', 0.30960,'DriftPass','EApertures',[Rx, Ry]);

% Straights added for achromat 11 scraper
str01925 = atdrift('str01925', 0.1925,'DriftPass','EApertures',[Rx, Ry]);
str03075 = atdrift('str03075', 0.3075,'DriftPass','EApertures',[Rx, Ry]);

%
% anti-bends scheme angle
%
q3   = 0.11; %deg

% Quadrupoles
% As some routines rely on the K value being set, equalize it to the
% PolynomB(2) value.  QF = 4.030928383463235
% QF      = atmultipole('QF'    , 0.150000, [0 0 0 0], [0  4.030076 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2);
% QFI     = atmultipole('QFI'   , 0.150000, [0 0 0 0], [-thU3/2 5.368969560153515 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2);

kQ1  =  3.981763113538452; %3.981763113511422; %  4.126576123705094;
kQ2  = -2.759093967069558; %-2.759093967201164; % -2.697078963828588;
kQ3  =  4.075364254221591; %4.075364254221591; %  4.210435344175177;
kQ4  =  4.367505587336566; %4.367505587336566; %  4.210435344175177;
kQ5  =  5.090518907635698; %5.090518907635698; % 5.469383459895300;
kQ6  =  5.937863598120707; %5.937863598120707; % 5.469383459895300;
kR1  =  4.504790935817660; %4.504790935817660; % 5.387917968472925;
kR2  =  6.252924524377427; %6.252924524377427; % 5.387917968472925;
Q1   = atmultipole('Q1' , 0.250000, [0 0 0 0], [0 kQ1 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[Rx, Ry]);  
Q2 	 = atmultipole('Q2' , 0.250000, [0 0 0 0], [0 kQ2 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[Rx, Ry]);   
Q3   = atmultipole('Q3' , 0.150000, [0 0 0 0], [0 kQ3 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[Rx, Ry]);   
Q4   = atmultipole('Q4' , 0.150000, [0 0 0 0], [0 kQ4 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[Rx, Ry]);
Q5   = atmultipole('Q5' , 0.150000, [0 0 0 0], [0 kQ5 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[Rx, Ry]);
Q6   = atmultipole('Q6' , 0.150000, [0 0 0 0], [0 kQ6 0 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[Rx, Ry]);
R1   = atsbend('R1' , 0.150000, -2*q3*pi/180, kR1, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
R2   = atsbend('R2' , 0.150000, -2*q3*pi/180, kR2, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
R1.K = R1.PolynomB(2);
R2.K = R2.PolynomB(2);
Q4.K = Q4.PolynomB(2);

Q3.K = Q3.PolynomB(2);
Q1.K = Q1.PolynomB(2);
Q2.K = Q2.PolynomB(2);
Q5.K = Q5.PolynomB(2);
Q6.K = Q6.PolynomB(2);

% -----
% SEXTs
% -----
k2S1 = -1.812366631091285e+02; %-1.806003401613026e+02; %-1.654986782091373e+02;
k2S2 =  3.079258122273534e+02; %3.079135487794222e+02; %2.934700000000000e+02;
k2S3 = -2.364877920365561e+02; %-2.364877920365561e+02; %-2.324930372170136e+02;
k2S4 =  2.765371984786756e+02; %2.695493861890224e+02;
k2S5 =  2.314201997403034e+02; %2.276491262143600e+02;

S1    = atmultipole('S1'  , 0.10000, [0 0 0 0], [0 0 k2S1 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[Rx, Ry]);   
S2h   = atmultipole('S2'  , 0.10000/2, [0 0 0 0], [0 0 k2S2 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[Rx, Ry]);   
S2    = cellcat(S2h,S2h);
%atmultipole('S2'  , 0.10000, [0 0 0 0], [0 0 k2S2 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[Rx, Ry]);   
S3    = atmultipole('S3'  , 0.10000, [0 0 0 0], [0 0 k2S3 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[Rx, Ry]);
S4h   = atmultipole('S4'  , 0.10000/2, [0 0 0 0], [0 0 k2S4 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[Rx, Ry]);  
S4    = cellcat(S4h,S4h); 
%atmultipole('S4'  , 0.10000, [0 0 0 0], [0 0 k2S4 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[Rx, Ry]);  
S5h   = atmultipole('S5'  , 0.10000/2, [0 0 0 0], [0 0 k2S5 0], 'StrMPoleSymplectic4Pass','MaxOrder',2,'EApertures',[Rx, Ry]); 
S5    = cellcat(S5h,S5h);

% ----
% OCTs
% ----
O1  = atmultipole('O1' , 0.10000, [0 0 0 0], [0 0 0 -1.801924271336929e+03], 'StrMPoleSymplectic4Pass','MaxOrder',3,'EApertures',[Rx, Ry]);  
O2  = atmultipole('O2' , 0.10000, [0 0 0 0], [0 0 0  1.981650239064280e+03], 'StrMPoleSymplectic4Pass','MaxOrder',3,'EApertures',[Rx, Ry]); 
O3  = atmultipole('O3' , 0.10000, [0 0 0 0], [0 0 0  5.737018706282323e+02], 'StrMPoleSymplectic4Pass','MaxOrder',3,'EApertures',[Rx, Ry]); 

%
% BENDS
% 
sU   = (1.094181+0.151199+0.151101+0.101861+0.001569+0.000089);

dip = -1.055010351848936; %-1.041250302267119;
seg = 5;
D0_U1d  = atsbend('D2', 0.361890/seg, 1.094181*pi/180/seg, -0.864858/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
D0_U1  = cellcat(repmat(cellcat(D0_U1d),1,seg));
DF1_U1 = atsbend('D2', 0.050000, 0.151199*pi/180, -0.864908/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DF2_U1 = atsbend('D2', 0.050000, 0.151101*pi/180, -0.866059/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DF3_U1 = atsbend('D2', 0.050000, 0.101861*pi/180, -0.551829/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DF4_U1 = atsbend('D2', 0.050000, 0.001569*pi/180,  0.011759/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DF5_U1 = atsbend('D2', 0.050000, 0.000089*pi/180, -0.000128/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);

D0_U2d   = atsbend('D2', 0.361890/seg, 1.094181*pi/180/seg , -0.864858/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
D0_U2   = cellcat(repmat(cellcat(D0_U2d),1,seg));
DF1_U2 = atsbend('D2', 0.050000, 0.151199*pi/180 , -0.864908/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DF2_U2 = atsbend('D2', 0.050000, 0.151101*pi/180 , -0.866059/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DF3_U2 = atsbend('D2', 0.050000, 0.101861*pi/180 , -0.551829/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DF4_U2 = atsbend('D2', 0.050000, 0.001569*pi/180 ,  0.011759/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DF5_U2 = atsbend('D2', 0.050000, 0.000089*pi/180 , -0.000128/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);

D0_U3d = atsbend('D3', 0.361890/seg, 1.094181* (1+4*q3/sU) *pi/180/seg, -0.864858/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
D0_U3  = cellcat(repmat(cellcat(D0_U3d),1,seg));
DF1_U3 = atsbend('D3', 0.050000, 0.151199* (1+4*q3/sU) *pi/180, -0.864908/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DF2_U3 = atsbend('D3', 0.050000, 0.151101* (1+4*q3/sU) *pi/180, -0.866059/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DF3_U3 = atsbend('D3', 0.050000, 0.101861* (1+4*q3/sU) *pi/180, -0.551829/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DF4_U3 = atsbend('D3', 0.050000, 0.001569* (1+4*q3/sU) *pi/180,  0.011759/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DF5_U3 = atsbend('D3', 0.050000, 0.000089* (1+4*q3/sU) *pi/180, -0.000128/-0.864858 * dip, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);

dipm =  -0.727588415518268; %-0.737696955554395;
DS6 = atsbend('D1', 0.050000, 0.001070*pi/180,  0.006608/ -0.870701 * dipm, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DS5 = atsbend('D1', 0.050000, 0.050729*pi/180, -0.271428/ -0.870701 * dipm, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DS4 = atsbend('D1', 0.050000, 0.074672*pi/180, -0.425119/ -0.870701 * dipm, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DS3 = atsbend('D1', 0.050000, 0.076248*pi/180, -0.426048/ -0.870701 * dipm, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DS2 = atsbend('D1', 0.050000, 0.114983*pi/180, -0.584884/ -0.870701 * dipm, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DS1 = atsbend('D1', 0.050000, 0.152049*pi/180, -0.870351/ -0.870701 * dipm, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DS0d = atsbend('D1', 0.204240/seg, 0.621695*pi/180/seg, -0.870701/ -0.870701 * dipm, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DS0 = cellcat(repmat(cellcat(DS0d),1,seg));


DM1 = atsbend('D1', 0.050000, 0.152220*pi/180, -0.870751/ -0.870701 * dipm, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DM2 = atsbend('D1', 0.050000, 0.152122*pi/180, -0.871910/ -0.870701 * dipm, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DM3 = atsbend('D1', 0.050000, 0.102549*pi/180, -0.555557/ -0.870701 * dipm, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DM4 = atsbend('D1', 0.050000, 0.001579*pi/180,  0.011839/ -0.870701 * dipm, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DM5 = atsbend('D1', 0.050000, 0.000090*pi/180, -0.000129/ -0.870701 * dipm, 'BndMPoleSymplectic4Pass','EApertures',[Rx, Ry]);
DIP1  = cellcat(DF5_U1, DF4_U1, DF3_U1, DF2_U1, DF1_U1, D0_U1);
DIP2  = cellcat(DF5_U2, DF4_U2, DF3_U2, DF2_U2, DF1_U2, D0_U2);
DIP3  = cellcat(DF5_U3, DF4_U3, DF3_U3, DF2_U3, DF1_U3, D0_U3);

DIPm = cellcat(DS6, DS5, DS4, DS3, DS2, DS1, DS0, DM1, DM2, DM3, DM4, DM5);

ShortStr_RF     = cellcat( str0100, corr_ac, str0362, CAV, str0462 );           % Preliminary position for fast corrector, mounted on taper

% BPM     = atmarker('BPM','IdentityPass');
BPM     = atmonitor('BPM', 'IdentityPass','EApertures',[Rx, Ry]);
BPMbloc = cellcat( str0025, BPM, str0025);

CMh     = atcorrector('CMh',0.0,[ 0 0 ],'CorrectorPass','iscorH','H');
CMv     = atcorrector('CMv',0.0,[ 0 0 ],'CorrectorPass','iscorV','V');
% CMh     = atmultipole('CMh',0.0,[0 0],[0 0],'StrMPoleSymplectic4Pass','MaxOrder',2); 
% CMv     = atmultipole('CMv',0.0,[0 0],[0 0],'StrMPoleSymplectic4Pass','MaxOrder',2); 

CORR    = cellcat( str0025, CMv, str0050, CMh, str0025);
CORRS   = cellcat( str0075, CMh, str0025);


PSM     = cellcat( str0150, str0150);
KI      = cellcat( str0150, str0150);
PV      = cellcat( str0075, str0075);
SEPT    = cellcat( stra500, stra500, stra025);

SQFM = cellcat(Q3, str0075, S2, strx013, BPMbloc, strx013, Q4, str0100, CORR);
SQFO = cellcat(Q5, str0075, S4, strx013, BPMbloc, strx013, Q6, str0100, CORR);
SQFI = cellcat(R1, str0075, S5, strx013, BPMbloc, strx013, R2, str0100, CORR);

%DIPUC = cellcat(SD, str0010, DIP, flip(DIP), str0010, SD);
DIPUC1 = cellcat(S3, str0010, DIP1, flip(DIP1), str0010, S3);
DIPUC2 = cellcat(S3, str0010, DIP2, flip(DIP2), str0010, S3);
DIPUC3 = cellcat(S3, str0010, DIP3, flip(DIP3), str0010, S3);
DIPUC4 = cellcat(S3, str0010, DIP2, flip(DIP2), str0010, S3);
DIPUC5 = cellcat(S3, str0010, DIP1, flip(DIP1), str0010, S3);

% OXX   = cellcat(str0050, OXXO, str0050);
% OXY   = cellcat(str0050, OXYO, str0050);
% OYY   = cellcat(str0050, OYYO, str0050);



LS      = cellcat( repmat(cellcat(str0500),1,4), str0321);
LS1A    = cellcat( str0500, str0377, SEPT, str0098, str0321);
LS2A    = cellcat( repmat(cellcat(str0500),1,3), str0252, PSM, str0269);
SS      = cellcat( repmat(cellcat(str0500),1,2), str0302);
SS1A    = cellcat( str0270, KI, str0150, str0280, str0302);
SS10B   = cellcat( str0420, str0075, PV, str0355, str0302);
UC1     = cellcat( SQFM, strx203, DIPUC1, strx403);
UC2     = cellcat( SQFO, strx203, DIPUC2, strx403);
UC3     = cellcat( SQFI, strx203, DIPUC3, strx203, flip(SQFI));
UC4     = cellcat( strx403, DIPUC4, strx203, flip(SQFO));
UC5     = cellcat( strx403, DIPUC5, strx203, flip(SQFM));
MC1     = cellcat( BPMbloc, str0058, flip(CORR), str0021, O1, str0025, Q1, str0025, O2,...
          str0100, Q2, strx006, DIPm, O3, strx093, CORRS, BPMbloc, str0020,...
          S1);
MC2     = cellcat( BPMbloc, str0058, flip(CORR), str0021, O1, str0025, Q1, str0025, O2,...
          str0100, Q2, strx006, DIPm, O3, strx093, CORR, BPMbloc, str0020,...
          S1);
%%% achrRF  = cellcat( LS, GS, MC1, GE, SS, UC1, UC2, UC3, UC4, UC5, ShortStr_RF, flip(MC2), flip(LS) );
achrRF  = cellcat( LS, GS,MC1,GE, SS, ...
                       GS,UC1,GE, ...
                       GS,UC2,GE, ...
                       GS,UC3,GE, ...
                       GS,UC4,GE, ...
                       GS,UC5,GE, ShortStr_RF, ...
                       GS,flip(MC2),GE, ...
                       flip(LS) );
achr    = cellcat( LS, GS,MC1,GE, SS, ...
                       GS,UC1,GE, ...
                       GS,UC2,GE, ...
                       GS,UC3,GE, ...
                       GS,UC4,GE, ...
                       GS,UC5,GE, SS, ...
                       GS,flip(MC2),GE, ...
                       flip(LS) );
%achr    = cellcat( LS, GS,MC1,GE, SS, ...
%                   UC1, UC2, UC3, UC4, UC5, flip(SS), flip(MC2), flip(LS));

%ring    = cellcat( achrRF, repmat(achr,1,19)); Periodicity=1; % , repmat(achr,1, 5));
 
ring    = cellcat( achrRF );Periodicity=20; 
%ring    = cellcat(achr); Periodicity=20;
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

% CAVITY CASE
ring = atSetRingProperties(ring, ...
    'name',LatticeFile, ...
    'Energy',E0, ...                                % AT2 energy unit is eV
    'HarmNumber',HarmNumber, ...
    'rf_frequency','nominal', ...
    'rf_voltage', V0*numel(findcells(ring,'FamName','cav')), ...
    'particle','relativistic', ...                  % AT tracking assumes beta = 1
    'periodicity',Periodicity,...
    'cavpts',findcells(ring,'FamName','cav')+1);    % +1 required as the information is added as the first element in the cell array

% NO CAVITY CASE
% % % ring = atSetRingProperties(ring, ...
% % %     'name',LatticeFile, ...
% % %     'Energy',E0, ...                                % AT2 energy unit is eV
% % %     'HarmNumber',HarmNumber, ...
% % %     'rf_frequency','nominal', ...
% % %     'rf_voltage', V0*numel(findcells(ring,'FamName','cav')), ...
% % %     'particle','relativistic', ...                  % AT tracking assumes beta = 1
% % %     'periodicity',Periodicity);



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
%evalin('caller','global THERING GLOBVAL ParamGroupDIP ParamGroupDIPM'); %%%removed MA 08032024
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
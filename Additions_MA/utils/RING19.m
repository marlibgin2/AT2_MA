function Acell = RING19
% emix = 33.0053 (pm)
%   U0 = 439.4937 (keV)
%   nu = 64.45017      24.249781
%   xi = -86.616954     -73.803901
%   stiffness = 3.327
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
CHAMBER =   aperture('CHAMBER', [-1, 1, -1, 1], 'EAperturePass');
FAMLIST{CHAMBER}.ElemData.Axes = [0.035 0.0055];                     %The given values are "half-gaps"

% RF definitions
CAV     = rfcavity('CAV', 0.378, V0, HarmNumber*PhysConst.c*RelBeta/L0, HarmNumber,'CavityPass');
landau  = drift('landau', 0.378, 'DriftPass');

%
% straights definition
%
delta = 0.12; % getting outer sextupoles S2 closer to [AB S1 AB] 
d1    = atdrift('d1', 0.054,'DriftPass'); % 0.15
d2    = atdrift('d2', 0.04 - delta,'DriftPass');% 0.175
d3    = atdrift('d3', 0.24 + delta,'DriftPass');%0.075
d4    = atdrift('d4', 0.556,'DriftPass');%0.01
d3s   = atdrift('d3s', 0.09+delta-0.15,'DriftPass');%0.075

dm1   = atdrift('dm1', 0.60,'DriftPass'); % 0.05
dm1a  = atdrift('dm1', 0.25,'DriftPass'); % 0.05
dm1b  = atdrift('dm1', 0.25,'DriftPass'); % 0.05

dm2   = atdrift('dm2', 0.22-0.1,'DriftPass'); % 0.02
dm3   = atdrift('dm3', 0.05,'DriftPass'); % 0.05
dstra = atdrift('dstra', 2.3032+7e-15+0.1,'DriftPass'); % 0.5

% delta = drift('d4', dlt1,'DriftPass');

%Fscal =1; % 12/20; 
% ---------
% Anti-Bend
% ---------
AB  =  atsbend('AB' , ...
              0.25, ...
             -0.0045, ... 
              4.062459314543996,...
              'BndMPoleSymplectic4Pass');  

ABm  =  atsbend('ABm' , ...
              0.25, ...
             -0.0045, ... 
              3.647424625563907,...
              'BndMPoleSymplectic4Pass');  

%tgb_th = 0.003932799272737 + ab_th; %0.31*pi/180  
TGB =  atsbend('TGB' , ...
              0.22, ...
              0.0095, ... 
              -3.642194102107883,...
              'BndMPoleSymplectic4Pass' ...   
            );  

TGBm =  atsbend('TGBm' , ...
              0.22, ...
              0.0095, ... 
             -3.690802884609243,...
              'BndMPoleSymplectic4Pass' ...   
            );  

%%%lgbh_th  = 0.985*pi/180; 
lgbh_th  = Xin(9)*pi/180; 
LGBh = atsbend('LGBh' , ...
              0.30, ...
              lgbh_th, ... 
              0, ...
              'BndMPoleSymplectic4Pass' ...   
            );

%lgbh_th =  0.039314788064924; %0.038222710618676;
%fac = [0.05 0.1 0.2 0.3 0.4];
fac = [0.01 0.02 0.2 0.2 0.5]*2.5/2.31; % *************
% % % fac = [f1 f2 f3 f4 f5]; % *************
%fac = [0.24  0.24  0.24  0.24  0.24 ];
%fac = [0.0098 0.02 0.2 0.2 0.5002]; % *************
% % % 
% % % LGB1 = sbend('LGB1' , ...
% % %               0.05, ...
% % %               lgbh_th * fac(1), ... 
% % %               0, ...5
% % %               0, ...
% % %               0, ...
% % %               'BndMPoleSymplectic4Pass' ...   
% % %             );
% % % LGB2 = sbend('LGB1' , ...
% % %               0.05, ...
% % %               lgbh_th * fac(2), ... 
% % %               0, ...
% % %               0, ...
% % %               0, ...
% % %               'BndMPoleSymplectic4Pass' ...   
% % %             );
% % % LGB3 = sbend('LGB1' , ...
% % %               0.05, ...
% % %               lgbh_th * fac(3), ... 
% % %               0, ...
% % %               0, ...
% % %               0, ...
% % %               'BndMPoleSymplectic4Pass' ...   
% % %             );
% % % LGB4 = sbend('LGB1' , ...
% % %               0.05, ...
% % %               lgbh_th * fac(4), ... 
% % %               0, ...
% % %               0, ...
% % %               0, ...
% % %               'BndMPoleSymplectic4Pass' ...   
% % %             );
% % % LGB5 = sbend('LGB1' , ...
% % %               0.05, ...
% % %               lgbh_th * fac(5), ... 
% % %               0, ...
% % %               0, ...
% % %               0, ...
% % %               'BndMPoleSymplectic4Pass' ...   
% % %             );
% % % 
%LGBh = [LGB1 LGB2 LGB3 LGB4 LGB5];
%TOTangle = sum(sum(lgbh_th.*fac) + tgb_th); 
ABangle  = -ab_th;
%disp(['CELL bending angle = ' num2str(TOTangle * 180/pi) ' (deg)'])
disp(['AB bending angle = ' num2str(ab_th * 180/pi) ' (deg)'])

%
% Sextupoles
%
% S1h   = multipole('S1h'   , 0.025000, [0 0 0 0], [0 0  1 0], 'StrMPoleSymplectic4Pass');
% S2    = multipole('S2'    , 0.050000, [0 0 0 0], [0 0 -1 0], 'StrMPoleSymplectic4Pass'); 

QM1  = quadrupole('QM1',0.10, 8.978490149765197,'StrMPoleSymplectic4Pass');
QM1a = quadrupole('QM1a',0.10, -0.213687321846813,'StrMPoleSymplectic4Pass');
QM2  = quadrupole('QM2',0.10, -4.782206124053459,'StrMPoleSymplectic4Pass');
QM3  = quadrupole('QM3',0.10, -2.383797160579114,'StrMPoleSymplectic4Pass');
QM4  = quadrupole('QM4',0.10, 8.595400266174163,'StrMPoleSymplectic4Pass');

S1h  = sextupole('S1h'   , 0.025000, +0, 'StrMPoleSymplectic4Pass');
S2   = sextupole('S2'    , 0.050000, -0, 'StrMPoleSymplectic4Pass');

TMEh     = [S1h, d1, AB, d2, S2, d3, TGB, d4, LGBh ];
TME      = [TMEh, reverse(TMEh) ];

matchC   = [d2,d2, d2, d2, AB, d2, d2, d2, d2, S1h];

%%%%% RING = TME; 
% pure TME
% RING     = [ reverse(LGBh) d4 TGB d3 S2 d2 AB d1 S1h repeat(TME,7)' S1h d1 AB d2 S2 d3 TGB d4 LGBh ]; 

% matching   
% RING     = [ QM1 dm1 QM2 dm2 reverse(LGBh) d4 TGB d3 S2 d2 AB d1 S1h ...
%              repeat(TME,7)' ...
%              S1h d1 AB d2 S2 d3 TGB d4 LGBh dm2 QM2 dm1 QM1]; 
%

% RING     = [ dstra QM1 dm1 QM2 dm2 reverse(LGBh) d4 TGBm dm3 QM3 d3s S2 d2 ABm d1 S1h ...
%                repeat(TME,7)' ...
%                S1h d1 ABm d2 S2 d3s QM3 dm3 TGBm d4 LGBh dm2 QM2 dm1 QM1 dstra]; 

RING     = [ dstra QM1 dm1a QM1a dm1b QM2 dm2 reverse(LGBh) d4 TGBm dm3 QM3 d3s S2 d2 ABm d1 S1h ...
             S1h d1 ABm d2 S2      d3s QM4 dm3      TGBm d4 LGBh LGBh d4 TGB d3 S2 d2 AB d1 S1h ... 
             repeat(TME,5)' ...
             S1h d1 AB d2 S2 d3 TGB d4 LGBh LGBh d4 TGBm      dm3 QM4 d3s     S2 d2 ABm d1 S1h ... 
             S1h d1 ABm d2 S2 d3s QM3 dm3 TGBm d4 LGBh dm2 QM2 dm1b QM1a dm1a  QM1 dstra]; 

%RING    = [ repeat(US,Nc)' ]

for i=1:length(FAMLIST)
    FAMLIST{i}.ElemData.Energy = Ekin;
    FAMLIST{i}.ElemData.NumIntSteps = 10; 
    FAMLIST{i}.ElemData.MaxOrder = 3;
end


%Build lattice
ELIST = RING;
buildlat(ELIST);

%Set fields 'Energy' and 'NumIntSteps'
% THERING = setcellstruct(THERING,'Energy',1:length(THERING), GLOBVAL.E0);
% for i=1:length(THERING), THERING{i}.NumIntSteps = 20; end

% Set up global parameters in the calling workspace, then exit
Scell = findspos(THERING,length(THERING)+1);
[nu, ksi] = tunechrom(THERING,1e-4,[18.25 2.2],'chrom',1e-8);
BBB = findcells(THERING,'Class','Bend');
Acell = 0;
for i=1:length(BBB)
  Acell = Acell + THERING{BBB(i)}.BendingAngle;
end
disp(['    cell length:   ' num2str(Scell,'%11.6g') ' m']);
disp(['     cell angle:   '  num2str(Acell,'%11.6g') ' (rad) / ' num2str(Acell*180/pi,'%11.6g') ' (deg)']);
disp(['  20xcell angle:   '  num2str(20*Acell,'%11.6g') ' (rad) / ' num2str(20*Acell*180/pi,'%11.6g') ' (deg)']);
%disp(['    R3 betatron tunes: ' num2str(nu(1),'%3.6g') ' / ' num2str(nu(2),'%3.6g') '  (H/V)']);
%disp(['    R3 chromaticities: ' num2str(ksi(1),'%3.6g') ' / ' num2str(ksi(2),'%3.6g') '  (H/V)']);
evalin('caller','global THERING FAMLIST GLOBVAL ParamGroupDIP ParamGroupDIPM');
disp('    Loading completed...');


function Acell = mTME19_match(Xin)
%  mTME9(4.96,.826); plotbetaeta; ats = atsummary; ats.naturalEmittance
%
%  X0 =  3.60   -3.50    0.0059    0.0070    0.150    0.175    0.075    0.010 // emix = 75.32
%  X1 =  4.00   -4.00    0.0040    0.0035    0.100    0.175    0.075    0.020 // emix = 26.05  
%  X2 = 4.0000   -4.00    0.0037    0.0031    0.100   0.150    0.070
%  0.030    0.9474 // emis = 26.46 

clear GLOBAL FAMLIST THERING GLOBVALlgb
global FAMLIST THERING GLOBVAL  %NOTE: GLOBVAL will be made obsolete shortly.

load PC;        %Load physical constants

% disp(' ');
ab_k   = Xin(1);%
tgb_k  = Xin(2);%
ab_th  = Xin(3);%0.35
tgb_th = Xin(4);%0.49
dd1    = Xin(5);% 
dd2    = Xin(6);%
dd3    = Xin(7);
dd4    = Xin(8);  

abm_k  =  1.78164;
tgbm_k = -2.22222; 
if length(Xin)>9
abm_k  = Xin(10); %ab_k;
tgbm_k = Xin(11); %tgb_k; 
end

% dlt1   = Xin(5);
% dlt2   = Xin(6); 
% f1     = Xin(7);
% f2     = Xin(8);
% f3     = Xin(9);
% f4     = Xin(10);
% f5     = Xin(11);

GLOBVAL.E0 = 3000e6;
GLOBVAL.LatticeFile = 'mTME16';

RelGamma = GLOBVAL.E0/PC.me_eV;
RelBeta =  sqrt(1 - RelGamma^-2);
Ekin = GLOBVAL.E0 - PC.me_eV;
HarmNumber = 176;        % # of RF buckets in the ring
L0 = 528;
V0 = 1.8e6 / 6;

% Vacuum chamber definition
CHAMBER =   aperture('CHAMBER', [-1, 1, -1, 1], 'EAperturePass');
FAMLIST{CHAMBER}.ElemData.Axes = [0.035 0.0055];                     %The given values are "half-gaps"

% RF definitions
CAV     = rfcavity('CAV', 0.378, V0, HarmNumber*PC.c*RelBeta/L0, HarmNumber,'CavityPass');
landau  = drift('landau', 0.378, 'DriftPass');

% straights definition
% strd1    = drift('strd1', d1,'DriftPass');
% strd2    = drift('strd1', d2,'DriftPass');
% strdiff  = drift('strdiff', D-d1-d2,'DriftPass');
delta = 0.12; % getting outer sextupoles S2 closer to [AB S1 AB] 
d1    = drift('d1', dd1,'DriftPass'); % 0.15
d2    = drift('d2', dd2 - delta,'DriftPass');% 0.175
d3    = drift('d3', dd3 + delta,'DriftPass');%0.075
d4    = drift('d4', dd4,'DriftPass');%0.01
d3s   = drift('d3', dd3+delta-0.15,'DriftPass');%0.075

dm1   = drift('dm1', 0.60,'DriftPass'); % 0.05
dm1a   = drift('dm1', 0.25,'DriftPass'); % 0.05
dm1b   = drift('dm1', 0.25,'DriftPass'); % 0.05

dm2   = drift('dm2', 0.22-0.1,'DriftPass'); % 0.02
dm3   = drift('dm3', 0.05,'DriftPass'); % 0.05
dstra = drift('dstra', 2.3032+7e-15+0.1,'DriftPass'); % 0.5

% delta = drift('d4', dlt1,'DriftPass');

%Fscal =1; % 12/20; 
% ---------
% Anti-Bend
% ---------
%ab_th =0.009558086209037; %.85* 0.78*pi/180 * Fscal;
AB  =  sbend('AB' , ...
              0.25, ...
             -ab_th, ... 
              0, ...
              0, ...
             +ab_k,...
              'BndMPoleSymplectic4Pass');  

ABm  =  sbend('ABm' , ...
              0.25, ...
             -ab_th, ... 
              0, ...
              0, ...
             +abm_k,...
              'BndMPoleSymplectic4Pass');  

%tgb_th = 0.003932799272737 + ab_th; %0.31*pi/180  
TGB =  sbend('TGB' , ...
              0.22, ...
              tgb_th + ab_th, ... 
              0, ...
              0, ...
              tgb_k,...
              'BndMPoleSymplectic4Pass' ...   
            );  

TGBm =  sbend('TGBm' , ...
              0.22, ...
              tgb_th + ab_th, ... 
              0, ...
              0, ...
              tgbm_k,...
              'BndMPoleSymplectic4Pass' ...   
            );  

%%%lgbh_th  = 0.985*pi/180; 
lgbh_th  = Xin(9)*pi/180; 
LGBh = sbend('LGBh' , ...
              0.30, ...
              lgbh_th, ... 
              0, ...5
              0, ...
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

QM1  = quadrupole('QM1',0.10, -6.18973,'StrMPoleSymplectic4Pass');
QM1a = quadrupole('QM1a',0.10, .7,'StrMPoleSymplectic4Pass');
QM2  = quadrupole('QM2',0.10, 10.9241,'StrMPoleSymplectic4Pass');
QM3  = quadrupole('QM3',0.10, -5.58204,'StrMPoleSymplectic4Pass');
QM4  = quadrupole('QM4',0.10, 2.81122,'StrMPoleSymplectic4Pass');

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


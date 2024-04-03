function check_mag_specs(ring)
%
%
%
% global ring 
TGBi       = findcells(ring,'FamName','TGB');
TGB_angle  = ring{TGBi(1)}.BendingAngle;
TGB_length = ring{TGBi(1)}.Length;
TGB_K      = ring{TGBi(1)}.PolynomB(2);

rho_theta  = TGB_length;
inv_rho    = TGB_angle / rho_theta; 
TGB_B      = inv_rho/0.2998 * 3; 
TGB_g      = 3/0.2998 * TGB_K;
disp(['+ ----------------------------------------------------------'])
disp(['| TGB '])
disp(['| B     =  ' num2str(TGB_B) '  (T)  /  rho   =  ' num2str(1/inv_rho) ' (m)'])
disp(['| theta =  ' num2str(TGB_angle) '  (rad) / '  num2str(TGB_angle*180/pi) ' (deg)' ])
disp(['| K     =  ' num2str(TGB_K) '  (m-2)  /  B"/2 =  ' num2str(TGB_g) '  (T/m)'])
disp(['+ ----------------------------------------------------------'])

TGBmi       = findcells(ring,'FamName','TGBm');
TGBm_angle  = ring{TGBmi(1)}.BendingAngle;
TGBm_length = ring{TGBmi(1)}.Length;
TGBm_K      = ring{TGBmi(1)}.PolynomB(2);

rho_theta  = TGBm_length;
inv_rho    = TGBm_angle / rho_theta; 
TGBm_B      = inv_rho/0.2998 * 3; 
TGBm_g      = 3/0.2998 * TGBm_K;
disp(['+ ----------------------------------------------------------'])
disp(['| TGBm '])
disp(['| B     =  ' num2str(TGBm_B) '  (T)  /  rho   =  ' num2str(1/inv_rho) ' (m)'])
disp(['| theta =  ' num2str(TGBm_angle) '  (rad) / '  num2str(TGBm_angle*180/pi) ' (deg)' ])
disp(['| K     =  ' num2str(TGBm_K) '  (m-2)  /  B"/2 =  ' num2str(TGBm_g) '  (T/m)'])
disp(['+ ----------------------------------------------------------'])



LGBhi       = findcells(ring,'FamName','LGBh');
LGBh_angle  = ring{LGBhi(1)}.BendingAngle;
LGBh_length = ring{LGBhi(1)}.Length;
LGBh_K      = ring{LGBhi(1)}.PolynomB(2);

rho_theta  = LGBh_length;
inv_rho    = LGBh_angle / rho_theta; 
LGBh_B      = inv_rho/0.2998 * 3; 
LGBh_g      = 3/0.2998 * LGBh_K;
disp(['+ ----------------------------------------------------------'])
disp(['| LGBh '])
disp(['| B     =  ' num2str(LGBh_B) '  (T)  / rho   =  ' num2str(1/inv_rho) ' (m)'])
disp(['| theta =  ' num2str(LGBh_angle) '  (rad) / '  num2str(LGBh_angle*180/pi) ' (deg)' ])
disp(['| K     =  ' num2str(LGBh_K) '  (m-2)  /  B"/2 =  ' num2str(LGBh_g) '  (T/m)'])
disp(['+ ----------------------------------------------------------'])



ABi       = findcells(ring,'FamName','AB');
AB_angle  = ring{ABi(1)}.BendingAngle;
AB_length = ring{ABi(1)}.Length;
AB_K      = ring{ABi(1)}.PolynomB(2);

rho_theta  = AB_length;
inv_rho    = AB_angle / rho_theta; 
AB_B      = inv_rho/0.2998 * 3; 
AB_g  = 3/0.2998 * AB_K;
disp(['+ ----------------------------------------------------------'])
disp(['| AB '])
disp(['| B     =  ' num2str(AB_B) '  (T)  / rho   =  ' num2str(1/inv_rho) ' (m)'])
disp(['| theta =  ' num2str(AB_angle) '  (rad) / '  num2str(AB_angle*180/pi) ' (deg)' ])
disp(['| K     =  ' num2str(AB_K) '  (m-2)  /  B"/2 =  ' num2str(AB_g) '  (T/m)'])
disp(['+ ----------------------------------------------------------'])


ABmi       = findcells(ring,'FamName','ABm');
ABm_angle  = ring{ABmi(1)}.BendingAngle;
ABm_length = ring{ABmi(1)}.Length;
ABm_K      = ring{ABmi(1)}.PolynomB(2);

rho_theta  = ABm_length;
inv_rho    = ABm_angle / rho_theta; 
ABm_B      = inv_rho/0.2998 * 3; 
ABm_g      = 3/0.2998 * ABm_K;
disp(['+ ----------------------------------------------------------'])
disp(['| ABm '])
disp(['| Bm    =  ' num2str(ABm_B) '  (T)  / rho   =  ' num2str(1/inv_rho) ' (m)'])
disp(['| theta =  ' num2str(ABm_angle) '  (rad) / '  num2str(ABm_angle*180/pi) ' (deg)' ])
disp(['| K     =  ' num2str(ABm_K) '  (m-2)  /  B"/2 =  ' num2str(ABm_g) '  (T/m)'])
disp(['+ ----------------------------------------------------------'])

%
% matching quads
%
QM1i      = findcells(ring,'FamName','QM1');
QM1_length = ring{QM1i(1)}.Length;
QM1_K      = ring{QM1i(1)}.PolynomB(2);

QM1_g  = 3/0.2998 * QM1_K;
disp(['+ ----------------------------------------------------------'])
disp(['| QM1 - L = ' num2str(QM1_length) ' (m)' ])
disp(['| K     =  ' num2str(QM1_K) '  (m-2)  /  B"/2 =  ' num2str(QM1_g) '  (T/m)'])
disp(['+ ----------------------------------------------------------'])

QM1ai      = findcells(ring,'FamName','QM1a');
QM1a_length = ring{QM1ai(1)}.Length;
QM1a_K      = ring{QM1ai(1)}.PolynomB(2);

QM1a_g      = 3/0.2998 * QM1a_K;
disp(['+ ----------------------------------------------------------'])
disp(['| QM1a - L = ' num2str(QM1a_length) ' (m)' ])
disp(['| K     =  ' num2str(QM1a_K) '  (m-2)  /  B"/2 =  ' num2str(QM1a_g) '  (T/m)'])
disp(['+ ----------------------------------------------------------'])


QM2i      = findcells(ring,'FamName','QM2');
QM2_length = ring{QM2i(1)}.Length;
QM2_K      = ring{QM2i(1)}.PolynomB(2);

QM2_g  = 3/0.2998 * QM2_K;
disp(['+ ----------------------------------------------------------'])
disp(['| QM2 - L = ' num2str(QM2_length) ' (m)'])
disp(['| K     =  ' num2str(QM2_K) '  (m-2)  /  B"/2 =  ' num2str(QM2_g) '  (T/m)'])
disp(['+ ----------------------------------------------------------'])

QM3i      = findcells(ring,'FamName','QM3');
QM3_length = ring{QM3i(1)}.Length;
QM3_K      = ring{QM3i(1)}.PolynomB(2);

QM3_g      = 3/0.2998 * QM3_K;
disp(['+ ----------------------------------------------------------'])
disp(['| QM3 - L = ' num2str(QM3_length) ' (m)'])
disp(['| K     =  ' num2str(QM3_K) '  (m-2)  /  B"/2 =  ' num2str(QM3_g) '  (T/m)'])
disp(['+ ----------------------------------------------------------'])

QM4i      = findcells(ring,'FamName','QM4');
QM4_length = ring{QM4i(1)}.Length;
QM4_K      = ring{QM4i(1)}.PolynomB(2);

QM4_g      = 3/0.2998 * QM4_K;
disp(['+ ----------------------------------------------------------'])
disp(['| QM4 - L = ' num2str(QM4_length) ' (m)'])
disp(['| K     =  ' num2str(QM4_K) '  (m-2)  /  B"/2 =  ' num2str(QM4_g) '  (T/m)'])
disp(['+ ----------------------------------------------------------'])


S1hi       = findcells(ring,'FamName','S1h');
S1h_length = ring{S1hi(1)}.Length;
S1h_K2     = ring{S1hi(1)}.PolynomB(3);

S1h_g2     = 3/0.2998 * S1h_K2;
disp(['+ ----------------------------------------------------------'])
disp(['| S1h - L = ' num2str(S1h_length) ' (m)'])
disp(['| K2      =  ' num2str(S1h_K2) '  (m-3)  /  B"/2 =  ' num2str(S1h_g2) '  (T/m2)'])
disp(['+ ----------------------------------------------------------'])


S2i       = findcells(ring,'FamName','S2');
S2_length = ring{S2i(1)}.Length;
S2_K2     = ring{S2i(1)}.PolynomB(3);

S2_g2     = 3/0.2998 * S2_K2;
disp(['+ ----------------------------------------------------------'])
disp(['| S2 - L = ' num2str(S2_length) ' (m)'])
disp(['| K2      =  ' num2str(S2_K2) '  (m-3)  /  B"/2 =  ' num2str(S2_g2) '  (T/m2)'])
disp(['+ ----------------------------------------------------------'])

HSi       = findcells(ring,'FamName','HS');
HS_length = ring{HSi(1)}.Length;
HS_K2     = ring{HSi(1)}.PolynomB(3);

HS_g2     = 3/0.2998 * HS_K2;
disp(['+ ----------------------------------------------------------'])
disp(['| HS - L = ' num2str(HS_length) ' (m)'])
disp(['| K2      =  ' num2str(HS_K2) '  (m-3)  /  B"/2 =  ' num2str(HS_g2) '  (T/m2)'])
disp(['+ ----------------------------------------------------------'])

HSFi       = findcells(ring,'FamName','HSF');
HSF_length = ring{HSFi(1)}.Length;
HSF_K2     = ring{HSFi(1)}.PolynomB(3);

HSF_g2     = 3/0.2998 * HSF_K2;
disp(['+ ----------------------------------------------------------'])
disp(['| HSF - L = ' num2str(HSF_length) ' (m)'])
disp(['| K2      =  ' num2str(HSF_K2) '  (m-3)  /  B"/2 =  ' num2str(HSF_g2) '  (T/m2)'])
disp(['+ ----------------------------------------------------------'])

HSDi       = findcells(ring,'FamName','HSD');
HSD_length = ring{HSDi(1)}.Length;
HSD_K2     = ring{HSDi(1)}.PolynomB(3);

HSD_g2     = 3/0.2998 * HSD_K2;
disp(['+ ----------------------------------------------------------'])
disp(['| HSD - L = ' num2str(HSD_length) ' (m)'])
disp(['| K2      =  ' num2str(HSD_K2) '  (m-3)  /  B"/2 =  ' num2str(HSD_g2) '  (T/m2)'])
disp(['+ ----------------------------------------------------------'])

OXXOi       = findcells(ring,'FamName','OXXO');
OXXO_length = ring{OXXOi(1)}.Length;
OXXO_K3     = ring{OXXOi(1)}.PolynomB(4);

OXXO_g3     = 3/0.2998 * OXXO_K3;
disp(['+ ----------------------------------------------------------'])
disp(['| OXXO - L = ' num2str(OXXO_length) ' (m)'])
disp(['| K2      =  ' num2str(OXXO_K3) '  (m-4) ']) %/  B"/2 =  ' num2str(HSD_g2) '  (T/m2)
disp(['+ ----------------------------------------------------------'])


end


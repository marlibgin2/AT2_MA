abmi  = findcells(RING_matched_optconstr,'FamName','ABm');
abi   = findcells(RING_matched_optconstr,'FamName','AB');
tgbi  = findcells(RING_matched_optconstr,'FamName','TGB');
tgbmi = findcells(RING_matched_optconstr,'FamName','TGBm');

qm0i  = findcells(RING_matched_optconstr,'FamName','QM0');
qm1i  = findcells(RING_matched_optconstr,'FamName','QM1');
qm1ai = findcells(RING_matched_optconstr,'FamName','QM1a');
qm2i  = findcells(RING_matched_optconstr,'FamName','QM2');
qm3i  = findcells(RING_matched_optconstr,'FamName','QM3');
qm4i  = findcells(RING_matched_optconstr,'FamName','QM4');

tgbm_A = RING_matched_optconstr{tgbmi(1)}.BendingAngle;
tgb_A  = RING_matched_optconstr{tgbi(1)}.BendingAngle;
abm_A  = RING_matched_optconstr{abmi(1)}.BendingAngle;
ab_A   = RING_matched_optconstr{abi(1)}.BendingAngle;

tgbm_K = RING_matched_optconstr{tgbmi(1)}.PolynomB(2);
tgb_K  = RING_matched_optconstr{tgbi(1)}.PolynomB(2);
abm_K  = RING_matched_optconstr{abmi(1)}.PolynomB(2);
ab_K   = RING_matched_optconstr{abi(1)}.PolynomB(2);

qm0_K  = RING_matched_optconstr{qm0i(1)}.PolynomB(2);
qm1_K  = RING_matched_optconstr{qm1i(1)}.PolynomB(2);
qm1a_K = RING_matched_optconstr{qm1ai(1)}.PolynomB(2);
qm2_K  = RING_matched_optconstr{qm2i(1)}.PolynomB(2);
qm3_K  = RING_matched_optconstr{qm3i(1)}.PolynomB(2);
qm4_K  = RING_matched_optconstr{qm4i(1)}.PolynomB(2);


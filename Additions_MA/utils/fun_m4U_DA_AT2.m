function penalty = fun_m4U_DA_AT2(PAR)


%
% 1st match the new tunw
%
NUtgt = [PAR(1) PAR(2)];
save('TUNEPOINT.mat','NUtgt');
[RING, RING_matched_optconstr, penalty, dmin] = runtotest_atmatch_m4U_b1_3_1__AT2;
RING_matched_optconstr__fitC = atfitchrom(RING_matched_optconstr, [1,1]/20,'S1','S3');

%
% alter the octupoles
%
oxxoi = findcells(RING,'FamName','OXXO');
oxyoi = findcells(RING,'FamName','OXYO');
oyyoi = findcells(RING,'FamName','OYYO');

for j = 1:2
        RING_matched_optconstr__fitC{oxxoi(j)}.PolynomB(4)=PAR(3);
        RING_matched_optconstr__fitC{oyyoi(j)}.PolynomB(4)=PAR(4);
        RING_matched_optconstr__fitC{oxyoi(j)}.PolynomB(4)=PAR(5);
end


DA = modelDA_sim_par(RING_matched_optconstr__fitC, 15e-3, 53, 2000, 0, 0.25e-3, 1.1);
[Area, ~] = calcDA_Area(DA);

penalty = -Area; 

end
clear all
nux = linspace(56.15,56.20,11);
nuy = linspace(16.11,16.36,11);
for i = 1:11
    for j = 1:11
        NUtgt = [nux(i), nuy(j)];
        save('TUNEPOINT.mat','NUtgt')
        [RING, RING_matched_optconstr] = runtotest_atmatch_m4U_b1_3_1__AT2;
        RING_matched_optconstr__fitC = atfitchrom(RING_matched_optconstr, [1,1]/20,'S1','S3');
        Rout = RING_matched_optconstr__fitC; 
        generaXY_DAgrid;
        figure(333)
        saveas(gcf, ['DA_' num2str(i) '_' num2str(j) '.png'])
        saveas(gcf, ['DA_' num2str(i) '_' num2str(j) '.fig'])

        9
    end
end

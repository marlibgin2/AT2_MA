
%nux = linspace(55.180,55.200,5);
%nuy = linspace(15.455,15.475,5);

nux = linspace(55.05,55.45,11);
nuy = linspace(15.05,15.45,11);
for i = 1:5
    for j = 1:5
        NU = [nux(i), nuy(j)];
        RoutTNOOCT = atfittune(Rin,[nux(i), nuy(j)],'QFEND','QDEND','UseIntegerPart');        
        RoutTNOOCT = atfittune(RoutTNOOCT,[nux(i), nuy(j)],'QFEND','QDEND','UseIntegerPart');
        RoutTCNOOCT = atfitchrom(RoutTNOOCT, [1,1],'SFM','SDEND');
        [x0pos, y0pos, nuxpos, nuypos, diffuvec,AE] = DA_NU_FMA_par(RoutTCNOOCT,'fma.out',66*8*2,80,80,7e-3,7e-3); 
        figure(35)
        text(-6,6,['NU = ' num2str(nux(i),6) '/' num2str(nuy(j),6)])
        saveas(gcf,['TS_DA_' num2str(nux(i),6) '_' num2str(nuy(j),6) '.png'])
        saveas(gcf,['TS_DA_' num2str(nux(i),6) '_' num2str(nuy(j),6) '.fig'])
        figure(36)        
        text(0.02,0.45,['NU = ' num2str(nux(i),6) '/' num2str(nuy(j),6)])
        saveas(gcf,['TS_NU_' num2str(nux(i),6) '_' num2str(nuy(j),6) '.png'])
        saveas(gcf,['TS_NU_' num2str(nux(i),6) '_' num2str(nuy(j),6) '.fig'])
    end
end



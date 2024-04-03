function J = Chrom_vs_sext(Rin)

names   = {'SD','SDEND','SFM','SFO','SFI'};
sdi     = findcells(Rin,'FamName','SD');
sdendi  = findcells(Rin,'FamName','SDEND');
sfmi    = findcells(Rin,'FamName','SFM');
sfoi    = findcells(Rin,'FamName','SFO');
sfii    = findcells(Rin,'FamName','SFI');

[twiss, tune, chrom]  = twissring(Rin,0.0, 1:length(Rin)+1, 'chrom', 1e-7);
Xix0 = chrom(1); Xiy0 = chrom(2); 
disp(['baseline   -- Xix = '  num2str(Xix0) ' Xiy = ' num2str(Xiy0)])

%
% modify sdi
% 

VARi = {sdi; sdendi; sfmi; sfoi; sfii};
for j = 1:length(VARi)
for i=1:length(VARi{j})
    Ktmp(i,j) = Rin{VARi{j}(i)}.PolynomB(3) ;
    Rin{VARi{j}(i)}.PolynomB(3)   =     Ktmp(i,j) * (1.01) ;
end
[twiss, tune, chrom]  = twissring(Rin,0.0, 1:length(Rin)+1, 'chrom', 1e-7);
Xixp = chrom(1); Xiyp = chrom(2); 
disp([names{j} ' +1% ----- Xix = '  num2str(Xixp) ' Xiy = ' num2str(Xiyp)])
for i=1:length(VARi{j})
    Rin{VARi{j}(i)}.PolynomB(3)   =     Ktmp(i,j) * 0.99 ;
end
[twiss, tune, chrom]  = twissring(Rin,0.0, 1:length(Rin)+1, 'chrom', 1e-7);
Xixm = chrom(1); Xiym = chrom(2); 
disp([names{j} ' -1% ----- Xix = '  num2str(Xixm) ' Xiy = ' num2str(Xiym)])
Jx(j) = (Xixp-Xixm)/(Ktmp(i,j) * 0.02);
Jy(j) = (Xiyp-Xiym)/(Ktmp(i,j) * 0.02);
disp([names{j} '     -----  Jx = '  num2str(Jx) '    Jy = ' num2str(Jy)])
for i=1:length(VARi{j})
    Rin{VARi{j}(i)}.PolynomB(3) = Ktmp(i,j) ;
end

end

J.Jx = Jx; J.Jy = Jy; 

dXidK2    = [J.Jx; J.Jy];
invdXidK2 = pinv(dXidK2);

Jred.Jx = [Jx(1), Jx(5)];
Jred.Jy = [Jy(1), Jy(5)];

dXidK2red    = [Jred.Jx; Jred.Jy]; % reduced Jacobian
invdXidK2red = pinv(dXidK2red);    % operates only on the sexts used to correct chrom

% test
dXi = [-3; +2];
dK2 = invdXidK2 * dXi;

for j = 1:length(VARi)
for i=1:length(VARi{j})
    Rin{VARi{j}(i)}.PolynomB(3)   =   Rin{VARi{j}(i)}.PolynomB(3) + dK2(j) ;
end
end

[twiss, tune, chrom]  = twissring(Rin,0.0, 1:length(Rin)+1, 'chrom', 1e-7);
Xix_new = chrom(1); Xiy_new = chrom(2); 
disp([names{j} ' >>> ----- Xix_new = '  num2str(Xix_new) ' Xiy_new = ' num2str(Xiy_new)])

% -------------------
% restore std machine
% -------------------
for j=1:5
    for i=1:length(VARi{j})
        Rin{VARi{j}(i)}.PolynomB(3) = Ktmp(i,j) ;
    end
end
% test 2
% alter sdendi; sfmi; sfoi;
%
dKtest = [0; 0.5; 0.9; 1.2; 0];
for j=1:5
    for i=1:length(VARi{j})
        Rin{VARi{j}(i)}.PolynomB(3) = Ktmp(i,j) + dKtest(j) ;
    end
end

DXi = dXidK2 * dKtest;
[twiss, tune, chrom]  = twissring(Rin,0.0, 1:length(Rin)+1, 'chrom', 1e-7);
Xix_alt = chrom(1); Xiy_alt = chrom(2); 
disp([names{j} ' >>> ----- Xix_new = '  num2str(Xix_alt) ' Xiy_new = ' num2str(Xiy_alt)])
disp([names{j} ' >>> ----- dXix = '  num2str(DXi(1)) ' dXiy = ' num2str(DXi(2))])

dK2 = invdXidK2red * (-DXi)
dKtest = [dK2(1); 0; 0; 0; dK2(2)];
for j=1:5
    for i=1:length(VARi{j})
        Rin{VARi{j}(i)}.PolynomB(3) = Rin{VARi{j}(i)}.PolynomB(3) + dKtest(j) ;
    end
end

9
end

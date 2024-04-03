function penalty = fun_numax_AT2(X0)
global rrr
%rrr = numax_7BA_24_2_1;

hsi = findcells(rrr,'FamName','HS');
hsfi = findcells(rrr,'FamName','HSF');
hsdi = findcells(rrr,'FamName','HSD');

oxxoi  = findcells(rrr,'FamName','oxxo');
oxyoi  = findcells(rrr,'FamName','oxyo');
oyyoi  = findcells(rrr,'FamName','oyyo');

for i = 1:length(hsi)
    rrr{hsi(i)}.PolynomB(3) = X0(1);
end
for i = 1:length(hsfi)
    rrr{hsfi(i)}.PolynomB(3) = X0(2);
end
for i = 1:length(hsdi)
    rrr{hsdi(i)}.PolynomB(3) = X0(3);
end

for i = 1:length(oxxoi)
    rrr{oxxoi(i)}.PolynomB(4) = X0(4);
end
for i = 1:length(oxyoi)
    rrr{oxyoi(i)}.PolynomB(4) = X0(5);
end
for i = 1:length(oyyoi)
    rrr{oyyoi(i)}.PolynomB(4) = X0(6);
end

tic; [x_da, y_da, DA] = moga_DA(rrr,15,40,500,10e-3,3e-3,0.0);toc
%%% tic; [x_da, y_da, DA] = moga_DA(rrr,3,40,1000,7e-3,3e-3,0.0);toc % the triangular DA! 
disp(['DA = ' num2str(DA)])
figure(112); hold on; plot(x_da, y_da, '-o'); 
penalty = -DA; 

end%
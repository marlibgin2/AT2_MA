function [penalty, Rout] = fun_maxmed_OCT(X)

global Rin

oxxoi = findcells(Rin,'FamName','oxxo');
oxyoi = findcells(Rin,'FamName','oxyo');
oyyoi = findcells(Rin,'FamName','oyyo');

VARi = {oxxoi; oxyoi; oyyoi};

for j = 1:3
for i=1:length(VARi{j})
    Rin{VARi{j}(i)}.PolynomB(4) = X(j);
end
end

clear VARi;

nturn=66*8*2;
[x0pos, y0pos, nuxpos, nuypos, diffuvec,AE] = DA_NU_FMA_par(Rin,'fma.out',nturn,40,40,7e-3,7e-3)
penalty = AE; 
Rout=Rin; 

end
function R=VarySextFam(R,K2val,fam)
%
% functions returns a new ring with a different value for K1 of fam

indfam=findcells(R,'FamName',fam);
%*ones(size(indfam))
R=setcellstruct(R,'PolynomB',indfam,K2val*ones(size(indfam)),1,3); 

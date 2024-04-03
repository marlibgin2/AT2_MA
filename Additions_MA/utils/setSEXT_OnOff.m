function Rout = setSEXT_OnOff(Rin, onoff)
Rout =Rin; 
sdi     = findcells(Rin,'FamName','SD');
sdendi  = findcells(Rin,'FamName','SDEND');
sfmi    = findcells(Rin,'FamName','SFM');
sfoi    = findcells(Rin,'FamName','SFO');
sfii    = findcells(Rin,'FamName','SFI');

if nargin<2 
    onoff=0;
end
sextFAC = onoff; 
SX = [ -116.625229 170.000000  170.000000 174.000000 207.412038 ]*sextFAC;

for i = 1:length(sdi)
   Rout{sdi(i)}.PolynomB(3) = SX(1); 
end
for i = 1:length(sdendi)
   Rout{sdendi(i)}.PolynomB(3) = SX(2); 
end
for i = 1:length(sfmi)
   Rout{sfmi(i)}.PolynomB(3) = SX(3); 
end
for i = 1:length(sfoi)
   Rout{sfoi(i)}.PolynomB(3) = SX(4); 
end
for i = 1:length(sfii)
   Rout{sfii(i)}.PolynomB(3) = SX(5); 
end





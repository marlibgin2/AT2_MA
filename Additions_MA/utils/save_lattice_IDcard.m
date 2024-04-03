Lattice_Name = 'm4U_240114_b01_02_03_02__grd';


RING = m4U_240114_b01_02_03_02__grd;
period = 20;
ringpars = check_main_params(RING, period);
Description = 'B type lattice after a series of MOGA optimsation, tune scan and optimsation to find the tune-octupoles producing a large DA';
lattMode = 'b01';
PB = findcells(RING,'PolynomB');
All_fams={};
for i = 1:numel(PB)
    All_fams{i} = RING{PB(i)}.FamName;
end
All_fams = unique(All_fams');
ACHROMAT=RING; 
%save([Lattice_Name '.mat'], 'ringpars','Description','ACHROMAT','lattMode','Lattice_Name','All_fams')
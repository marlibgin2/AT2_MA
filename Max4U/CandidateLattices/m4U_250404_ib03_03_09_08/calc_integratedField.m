function [Tm T] = calc_integratedField(ACHRO,FamName,varargin)

BendAng = getoption(varargin,'BendAng',0); %rad
K1      = getoption(varargin,'K1',0); %
K2      = getoption(varargin,'K2',0); %
switch FamName    
    case 'S1'
        iS = findcells(ACHRO,'FamName','S1'); N=1;
    case 'S3'
        iS = findcells(ACHRO,'FamName','S3'); N=1;
    case 'S6'
        iS = findcells(ACHRO,'FamName','S6'); N=1;
end

Sk = 0; 
for i=1:N
    Sk = Sk + 10/2.998*3 * K1 * ACHRO{iS(i)}.Length; 
end
Tm = 0; T = Sk; 

end
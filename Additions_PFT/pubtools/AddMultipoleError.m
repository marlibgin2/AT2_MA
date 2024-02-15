function [ACHRO_err] = AddMultipoleError(ACHRO,an,bn,r0,famnames,n0,mscale)
% Adds a set of systematuc multipolar components to a given magnet fanily
%% inputs 
%    ACHROMAT: At 2.0 lattice to which errors ar to be added
%    an: row vetor  with normalized skew multipolar components
%    bn: row vector with normalized normal multipolar components
%    famnames: magnet families to which mutipoles will be added
%    n0: main multipole (n0=2 is quadrupole)
%    reference normalization radius 
%    scaling factor for all multipoles
%
nmults = max(size(an,2),size(bn,2));
bn(n0)=1.0E4/mscale;
an(n0)=0.0;

PolynomA    = zeros(1,nmults);
PolynomB    = zeros(1,nmults);
nfams = size(famnames,2);
ACHRO_err=ACHRO;
for i=1:nfams
    I_fam = find(atgetcells(ACHRO, 'FamName', famnames{i}));
    K_fam = atgetfieldvalues(ACHRO, I_fam, 'PolynomB', {n0});
    K0 = K_fam(1);  

    for n=1:nmults
        PolynomA(n) = an(n)*K0/(r0^(n-n0))*1E-4*mscale;
        PolynomB(n) = bn(n)*K0/(r0^(n-n0))*1E-4*mscale; 
    end
    nelem = size(I_fam,1);

    for j=1:nelem
       ACHRO_err{I_fam(j)}.PolynomA=PolynomA;
       ACHRO_err{I_fam(j)}.PolynomB=PolynomB;
       ACHRO_err{I_fam(j)}.MaxOrder = nmults-1;
    end
end

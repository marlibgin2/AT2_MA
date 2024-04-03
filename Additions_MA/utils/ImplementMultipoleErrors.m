
% ImplementMultipoleErrors 
RING = m4U_240114_b01_02_03_02__grd;

% -----------
% QUADRUPOLES
% -----------
Q1i = findcells(RING,'FamName','Q1');
Q2i = findcells(RING,'FamName','Q2');
Q3i = findcells(RING,'FamName','Q3');
Q4i = findcells(RING,'FamName','Q4');
Q5i = findcells(RING,'FamName','Q5');
Q6i = findcells(RING,'FamName','Q6');
Qi  = [Q1i Q2i Q3i Q4i Q5i Q6i];

an  = [0 0.12 0 0 0 -0.07 0 0 0 -0.03 0 0 0 -0.01]'*1e-4;
bn  = [0 1e4  0 0 0  0.14 0 0 0 -0.92 0 0 0 -0.12]'*1e-4; 

an0 = [0 0    0 0 0  0    0 0 0  0    0 0 0  0]'*1e-4;
bn0 = [0 1e4  0 0 0  0    0 0 0  0    0 0 0  0]'*1e-4;

Ean = an - an0;
Ebn = bn - bn0; 

rho0 = 7e-3; N=2; 
scf=1000;
[RINGm,Pnew,Panew]=AssignFieldErr(RING,Qi,N,rho0,scf*Ebn',scf*Ean');



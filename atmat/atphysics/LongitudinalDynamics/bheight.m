function bucketheight = bheight(V0,U0,alphac,E0,h)
% Calculates RF bucket height
%% Inputs
% Mandatory inputs
% V0: RF voltage [V]
% U0: energy loss per turn [eV]
% alphac : momentum compaction factor
% E0 : beam energy [eV]
% h : harmonic number
%
%% Outputs
% bucketheight: maximum relative momentum deviation 
%
% Usage examples
% E0=3E9; V0=1.0E6; U0=360E3; alphac=3E-4; h=176;
% bh=bheight(V0,U0,alphac,E0,h);

%% History
% PFT 2024/08/06: first version

q = V0/U0; % overvoltage
bucketheight = sqrt(2*U0/(pi*h*alphac*E0)*(sqrt(q^2-1)-acos(1/q)));


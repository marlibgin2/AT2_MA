function V0 = VvsBH(bh,U0,alphac,E0,h)
% Calculates the rf voltage needed to achieve a given
% rf bucket height
%
%% Mandatory Inputs
% bh: Desired Bucket Height
% U0: energy loss per turn [eV]
% alphac : momentum compaction factor
% E0 : beam energy [eV]
% h : harmonic number
%
%% Outputs
% V0 : RF voltage [V]
%
%% Usage examples
% bh=0.04;  U0=360E3; alphac=3E-4; E0=3E9;h=176;
% V0=VvsBH(bh,U0,alphac,E0,h);
% V0=VvsBH(0.04,360e3, 3e-4, 3e9,176);

%% History
% PFT 2024/08/06: function structure
% S Jena 2024/08/06 : solve for V0 for a desired bucket height using fzero

bucket_height = @(V0) bheight(V0, U0,alphac,  E0,h, bh);

V0_lb = 1.01* U0; 
V0_ub = 10 * U0;

options = optimset('TolX', 1e-9, 'Display', 'off');
V0 = fzero(bucket_height, [V0_lb, V0_ub], options);
end

function F = bheight(V0, U0,alphac,  E0,h, bh)
q = V0 / U0; % Overvoltage
bucketheight = sqrt(2*U0/(pi*h*alphac*E0)*(sqrt(q^2-1)-acos(1/q)));
F = bucketheight - bh;
end


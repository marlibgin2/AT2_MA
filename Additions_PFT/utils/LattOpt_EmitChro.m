function f = LattOpt_EmitChro(x,LatticeOptData)
% inputs
%       x: row vector with decision variables for Lattice Optimization.
%       LatticeOptData: structure with optimization configuration data, 
%       including the list of decision and objective variable and lattice.
%
% output: row vector of fitness functions to be minimized
%
%   

ACHRO    = LatticeOptData.ACHRO;
isdipole = LatticeOptData.isdipole;

ACHRO  = setDVs(2, ACHRO,LatticeOptData, x);

%
% Calculates Objective functions
%
try
    rpara = atsummary_fast(ACHRO,isdipole);
    Emitt = rpara.naturalEmittance*1E12;
    ChroX = rpara.chromaticity(1);
    ChroY = rpara.chromaticity(2);
    Jx    = rpara.damping(1);
    if ( (Emitt>0) && (Jx>0) && (Jx<3) )
       f(1) = Emitt;
       f(2) = sqrt(ChroX^2+ChroY^2);
     else
       f(1)=Inf;
       f(2)=Inf;
     end
catch ME
    fprintf('error in LattOpt_EmitChro:atsummary_fast for x = %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n', x);
    fprintf('Error message was:%s \n',ME.message);
    f(1)=Inf;
    f(2)=Inf;
end
end


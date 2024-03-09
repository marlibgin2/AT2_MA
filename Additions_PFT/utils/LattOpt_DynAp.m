function f = LattOpt_DynAp(x, LatticeOptData)
% inputx: 
%        x:row vector with decision variables for Lattice Optimization.
%        LatticeOptData: structure with optimization data, including the
%                        list of decision and objective variables
%
% output: row vector of objective functions to be minimized.
%
% The dynamic aperture area in real space (possibly scaled to 
% chosen beta functions is the objective function.
% Chromaticity is corrected to chosen target values with two
% chosen sextupole sextupole families. These are chosen when running
% the script m4 and are recorded in LattOptData.
%
%
% Parameters for dynamic aperture calculation
%
chroms0    = LatticeOptData.DAoptions.chroms0; % Target chromaticity for one superperiod
TolChrom   = LatticeOptData.DAoptions.TolChrom;% Chromaticity tolerances
Nitchro    = LatticeOptData.DAoptions.Nitchro; % Max n. iterations of chromaticty correction

chrom_fams = LatticeOptData.chrom_fams; % list of sextupole families to be used for chromaticity correction

ACHRO       = LatticeOptData.ACHRO;
isdipole    = LatticeOptData.isdipole;

ACHRO  = setDVs(2, ACHRO,LatticeOptData, x);
%
% Calculates Objective functions
%

%
% Verifies linear lattice parameters for optMode=Full
% 
try
    if (strcmpi(LatticeOptData.optMode,'Full'))
        rpara = atsummary_fast(ACHRO,isdipole);
        Emitt = rpara.naturalEmittance*1E12;
        etax  = rpara.etax;
        Jx    = rpara.damping(1);
        betax = rpara.beta0(1);
        betay = rpara.beta0(2);
        if ( (Emitt<=0) || (Jx<=0) && (Jx>=3) )
            f=0;
            return;
        end
    else
        etax = 0.0; %% temporary fix to maintain speed for "Nonlinear" optimization mode
        betax = NaN;
        betay = NaN;
    end
catch
    f=0;
    return
end
%
% Corrects chromaticity
%
try
   [ACHRO,~,~] = fitchroit(ACHRO, chrom_fams, chroms0, Nitchro, TolChrom);    
catch ME
   fprintf('Error in LattOpt_DynAp: chromaticity fit \n');
   fprintf('Error message was: %s \n',ME.message);
end
%
% Calculates dynamic aperture
%
try
    [DA,~]=calcDA_raw(ACHRO,LatticeOptData.DAoptions,etax,betax,betay);
    f=-DA;
catch ME
    fprintf('Error in LattOpt_DynAp: Dynamic Aperture calculation \n');
    fprintf('Error message was: %s \n',ME.message);    
end

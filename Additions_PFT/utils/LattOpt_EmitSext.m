function f = LattOpt_EmitSext(x,LatticeOptData)
%% Inputs
%       x: row vector with decision variables for Lattice Optimization.
%       LatticeOptData: structure with optimization configuration data, 
%       including the list of decision and objective variable and lattice.
%
%% Outputs: row vector of fitness functions to be minimized
%       [emittance sqrt(chrox^2 + chroy^2)]
%
%   
%%
ACHRO    = LatticeOptData.ACHRO;
isdipole = LatticeOptData.isdipole;

chroms0  = LatticeOptData.DAoptions.chroms0; % Target chromaticity for one superperiod
TolChrom = LatticeOptData.DAoptions.TolChrom;% Chromaticity tolerances
Nitchro  = LatticeOptData.DAoptions.Nitchro; % Max n. iterations of chromaticty correction
nsext_fams = LatticeOptData.nsext_fams;
chrom_fams = LatticeOptData.chrom_fams;
Isexts     = LatticeOptData.Isexts;
Iss = zeros(nsext_fams,1);
for i=1:nsext_fams
    Iss(i)=LatticeOptData.Isexts{i}(1);
end

ACHRO  = setDVs(2, ACHRO,LatticeOptData, x);

%
% Calculates Objective functions
%
try
    rpara = atsummary_fast(ACHRO,isdipole);
    Emitt = rpara.naturalEmittance*1E12;
    Jx    = rpara.damping(1);
    if ( (Emitt>0) && (Jx>0) && (Jx<3) )
       f(1) = Emitt;
       try
            [ACHRO,~,~] = fitchroit(ACHRO, chrom_fams, chroms0, Nitchro, TolChrom);    
            sexts=0;
            for i=1:nsext_fams
                sexts=sexts+(ACHRO{Iss(i)}.PolynomB(1,3))^2;
            end
            sexts=sqrt(sexts);
            f(2)=sexts;
       catch ME
            fprintf('Error in Latt_OptEmitSexts: chromaticity fit \n');
            fprintf('Error message was: %s \n',ME.message);
            f(2)=Inf;
       end
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


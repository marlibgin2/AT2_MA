function f = UCOpt_EmitSext(x,LatticeOptData)
% inputs
%       x: row vector with decision variables for Lattice Optimization.
%       LatticeOptData: structure with optimization configuration data, 
%       including the list of decision and objective variable and lattice.
%
% output: row vector of fitness functions to be minimized
%
%   
if (isfield(LatticeOptData,'Emit0'))
    Emit0=LatticeOptData.Emit0;
else
    Emit0=0;
end

UC       = LatticeOptData.UC;
isdipole = LatticeOptData.isdipoleUC;

chroms0  = LatticeOptData.DAoptions.chroms0; % Target chromaticity for one superperiod
TolChrom = LatticeOptData.DAoptions.TolChrom;% Chromaticity tolerances
Nitchro  = LatticeOptData.DAoptions.Nitchro; % Max n. iterations of chromaticty correction
nsext_fams = LatticeOptData.nsext_fams;
chrom_fams = LatticeOptData.chrom_fams;
Isexts     = LatticeOptData.IsextsUC;
Iss = zeros(nsext_fams,1);
for i=1:nsext_fams
    if (not(isempty(Isexts{i})))
        Iss(i)=Isexts{i}(1);
    else
        Iss(i)=[];
    end
end
Iss(Iss==0) = [];

UC  = setDVs(3, UC, LatticeOptData, x);

%
% Calculates Objective functions
%
try
    rpara = atsummary_fast(UC,isdipole);
    Emitt = rpara.naturalEmittance*1E12;
    ChroX = rpara.chromaticity(1);
    ChroY = rpara.chromaticity(2);
    Jx    = rpara.damping(1);
    if ( (Emitt>0) && (Jx>0) && (Jx<3) )
       f(1) = abs(Emitt-Emit0);
       try
            [UC,~,~] = fitchroit(UC, chrom_fams, [0 0], Nitchro, TolChrom);    
            sexts=0;
            for i=1:length(Iss)
                sexts=sexts+(UC{Iss(i)}.PolynomB(1,3))^2;
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
    fprintf('error in UCOpt_EmitChro:atsummary_fast for x = %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n', x);
    fprintf('Error message was:%s \n',ME.message);
    f(1)=Inf;
    f(2)=Inf;
end
end


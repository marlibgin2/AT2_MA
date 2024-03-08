function f = LattOpt_EmitRDT(x,LatticeOptData)
% inputs
%       x: row vector with decision variables for Lattice Optimization.
%       LatticeOptData: structure with optimization configuration data, 
%       including the list of decision and objective variable and lattice.
%
% output: row vector of fitness functions to be minimized: emittance and
% a penalty calculated from RDTs (following OPA) is
% returned. Chromaticity is corrected to (+1/+1) - for 20 achromats 
% with the sextupole families in chrom_fams. 
% 
%

% Parameters for RDT calculation
%
chroms0  = LatticeOptData.DAoptions.chroms0; % Target chromaticity for one superperiod
TolChrom = LatticeOptData.DAoptions.TolChrom;% Chromaticity tolerances
Nitchro  = LatticeOptData.DAoptions.Nitchro; % Max n. iterations of chromaticty correction
nsext_fams = LatticeOptData.nsext_fams;
Isexts     = LatticeOptData.Isexts;
Iss = zeros(nsext_fams,1);
for i=1:nsext_fams
    Iss(i)=LatticeOptData.Isexts{i}(1);
end

chrom_fams = LatticeOptData.chrom_fams; % list of sextupole families to be used for chromaticity correction

ACHRO       = LatticeOptData.ACHRO;
isdipole    = LatticeOptData.isdipole;

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
          k=0;   
          [~,chr]=tunechrom(ACHRO);
           while ( (abs(chr(1)-chroms0(1))>TolChrom(1))||...
                  (abs(chr(2)-chroms0(2))>TolChrom(2)) &&(k<Nitchro))
              ACHRO=atfitchrom(ACHRO,chroms0,chrom_fams{1},chrom_fams{2});
              [~,chr]=tunechrom(ACHRO);
              k=k+1;
           end
           try
                RDTs = computeRDT(ACHRO,1);
                f(2) = RDT_Penalty(RDTs);
                for i=1:nsext_fams
                    SS = ACHRO{Iss(i)}.PolynomB(3)^2;
                    f(2) = f(2) + SS*(0.1^2); % all sextupoles are 10 cm long
                end
           catch
                fprintf('Not possible to calculate RDTs\n');
                f(2)=Inf;       
            end
        catch
            disp('Chromaticity fitting not possible');
            f(2)=Inf;
       end
     else
       f(1)=Inf;
       f(2)=Inf;
     end
catch ME
    fprintf('error in at_summary for x = %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n', x);
    fprintf('Error message was:%s \n',ME.message);
    f(1)=Inf;
    f(2)=Inf;
end
end


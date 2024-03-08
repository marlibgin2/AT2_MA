function f = LattOpt_RDT(x,LatticeOptData)
% Calculates Penalty function based on Resonance Driving Terms (RDTs) 
%  
% inputs: 
%        x:row vector with decision variables 
%
%        LatticeOptData: structure with optimization data, including the
%                        list of decision and objective variable
%
% output: row vector of fitness functions to be minimized.
%
% A penalty function calculated from RDTs (following OPA) is
% returned. Chromaticity is corrected to (+1/+1) - for 20 achromats 
% with the sextupole families specified in chrom_fams. 
%
%
% Parameters for RDT calculation
%
chroms0    = LatticeOptData.DAoptions.chroms0; % Target chromaticity for one superperiod
TolChrom   = LatticeOptData.DAoptions.TolChrom;% Chromaticity tolerances
Nitchro    = LatticeOptData.DAoptions.Nitchro; % Max n. iterations of chromaticty correction
chrom_fams = LatticeOptData.chrom_fams; % list of sextupole families to be used for chromaticity correction
nsext_fams = LatticeOptData.nsext_fams;
Isexts     = LatticeOptData.Isexts;

ACHRO       = LatticeOptData.ACHRO;
isdipole    = LatticeOptData.isdipole;

ACHRO  = setDVs(2, ACHRO, LatticeOptData, x);
for i=1:nsext_fams
    Iss(i)=LatticeOptData.Isexts{i}(1);
end

%
% Calculates Objective functions
%
try
    rpara = atsummary_fast(ACHRO,isdipole);
    Emitt = rpara.naturalEmittance*1E12;
    Jx    = rpara.damping(1);
    if ( (Emitt>0) && (Jx>0) && (Jx<3) )
        try
          k=0;   
          [~,chr]=tunechrom(ACHRO);
           while ( (abs(chr(1)-chroms0(1))>TolChrom(1))||...
                  (abs(chr(2)-chroms0(2))>TolChrom(2)) &&(k<Nitchro))
              ACHRO=atfitchrom(ACHRO,chroms0,chrom_fams{1},chrom_fams{2});
              [~,chr]=tunechrom(ACHRO);
              k=k+1;
           end
        catch
            disp('Chromaticity fitting not possible');
       end
%
% Calculate RDTs
% 
       try
           RDTs = computeRDT(ACHRO,1);
           f    = RDT_Penalty(RDTs); 
           for i=1:nsext_fams
               SS = ACHRO{Iss(i)}.PolynomB(3)^2;
               f = f + SS*(0.1^2); % all sextupoles are 10 cm long
           end
       catch ME
           fprintf('Not possible to calculate RDTs\n');
           fprintf('Error message was:%s \n',ME.message);
           f=Inf;       
       end
    else
       f=Inf;
    end
catch ME
    fprintf('Error in atsummary_fast\n');
    fprintf('Error message was:%s \n',ME.message);
    f=Inf;
end

end



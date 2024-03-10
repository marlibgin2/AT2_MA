function f = RINGOpt_EmitDynAp(x, LatticeOptData)
%% input: 
%        x:row vector with decision variables for Lattice Optimization.
%        LatticeOptData: structure with optimization data, including the
%                        list of decision and objective variables
%
%        LatticeOptData.RINGGRD is the input lattice. This contains the
%        full ring, possibly with errors
%
%% output: row vector of fitness functions to be minimized.
% f(1) = Emittance
% f(2) = DA
%
% The dynamic aperture area in real space (possibly scaled to chosen 
% beta functions) and emittance are objective functions. 
% Chromaticity is corrected to chosen target values with two
% chose sextupole sextupole families. These are chosen when running
% the scripts max4_UpgradeStudies/m4U and are recorded in LattOptData.
%
%% Comments relevant for optmization configured with m4U.m
%     This can used for MOGA runs with optimization modes 'Linear' or 'Full'.
%     Optimization mode "Non-linear" intended for SOGA 
%
%% Parameters for dynamic aperture calculation
%
PC=load('PC.mat');      %to prevent matlab from complaining about variable name being the same as script name.
PhysConst = PC.PC;      %Load physical constants
%
ErrorModel = LatticeOptData.ErrorModel;
DAoptions  = LatticeOptData.DAoptions;
chroms0    = DAoptions.chroms0; % Target chromaticity for the whole ring
TolChrom   = DAoptions.TolChrom;% Chromaticity tolerances
Nitchro    = DAoptions.Nitchro; % Max n. iterations of chromaticty correction
TRmode     = DAoptions.TRmode; % tracking mode is 4d or 6d

chrom_fams = LatticeOptData.chrom_fams; % list of sextupole families to be used for chromaticity correction

ACHRO           = LatticeOptData.ACHRO;
isdipole        = LatticeOptData.isdipole;
RINGGRD         = LatticeOptData.RINGGRD;
%isdipoleRINGGRD = LatticeOptData.isdipoleRINGGRD;

ACHRO   = setDVs(2, ACHRO,LatticeOptData, x);

%% Calculates Objective functions
%
try
    rpara = atsummary_fast(ACHRO,isdipole);
    Emitt = rpara.naturalEmittance*1E12;
    etax  = rpara.etax;
    Jx    = rpara.damping(1);
    betax = rpara.beta0(1);
    betay = rpara.beta0(2);
    if ( (Emitt>0) && (Jx>0) && (Jx<3) )
       f(1) = Emitt;
       try
            [ACHRO,~,~] = fitchroit(ACHRO, chrom_fams, chroms0, Nitchro, TolChrom);    
       catch ME
            fprintf('Error in Latt_OptEmitDynAp: chromaticity fit \n');
            fprintf('Error message was: %s \n',ME.message);
       end
%
% Calculates dynamic aperture
%
       XAll=getAllfams(2,ACHRO,LatticeOptData);
       RINGGRD = setAllfams(6,RINGGRD,LatticeOptData,XAll);
       if (isstruct(ErrorModel))
        RINGGRD = applyErrorModel(RINGGRD,ErrorModel);
        RINGGRD = calcOrb(RINGGRD,'correct');
       end
       if(strcmpi(TRmode,'4d'))
           RINGGRD=atdisable_6d(RINGGRD);
       else
            ats=atsummary(RINGGRD);
            DAoptions.z0 = PhysConst.c*(ats.syncphase-pi)/(2*pi*ats.revFreq*ats.harmon); % choose the synchronous phase
       end
       try
            [DA,~]=calcDA_raw(RINGGRD,DAoptions,etax,betax,betay);
            f(2)=-DA;
       catch ME
            fprintf('Error in LattOpt_EmitDynAp: Dynamic Aperture calculation \n');
            fprintf('Error message was: %s \n',ME.message);    
       end
    else
       f(1)=Inf;
       f(2)=0;
    end
catch
    f(1)=Inf;
    f(2)=0;
end
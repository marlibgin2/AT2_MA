function f = LattOpt_EmitDynAp(x, LatticeOptData)
%% input: 
%        x:row vector with decision variables for Lattice Optimization.
%        LatticeOptData: structure with optimization data, including the
%                        list of decision and objective variables
%
%        LatticeOptData.ACHRO is the input lattice: this is one achromat.
%
%% output: row vector of fitness functions to be minimized.
% f(1) = Emittance
% f(2) = DA
%
%
% The dynamic aperture area of in real space (possibly scaled to chosen 
% beta functions) and emittance are objective functions. 
% Chromaticity is corrected to chosen target values with two
% chose sextupole sextupole families. These are chosen when running
% the scripts max4_UpgradeStudies/m4U and are recorded in LattOptData.
%
%% Comments relevant for optimization configured with max4_UpgradeStudies.m (Legacy)
% For revBmode = "All", this function may be used in two optmization modes 
%     SEXT : only the three sextupole family strengths above are
%            variables: DVs DV1 to DV3; (this is actually more appropriate
%            for a SOGA run - since the emittance is not changed by the 
%            sextupole settings)
%
%     CHRO : 7 decision variables define the linear lattice DV01 too DV07 , 
%            whereas DV08 to 10 are the three sextupoles families above.
%
% For revBmode="U3", this function may be used for one
%               optimization mode only (at the moment) :
%
%     CHRO: 14 decision variables, 11 defining the linear lattice and the
%     three sextupole families above.
%
%% Comments relevant for optmization configured with m4U.m
%     This can used for MOGA runs with optimization modes 'Linear' or 'Full'.
%     Optimization mode "Non-linear" intended for SOGA 
%
%% Parameters for dynamic aperture calculation
%
DAoptions  = LatticeOptData.DAoptions;
chroms0    = DAoptions.chroms0; % Target chromaticity for one superperiod
TolChrom   = DAoptions.TolChrom;% Chromaticity tolerances
Nitchro    = DAoptions.Nitchro; % Max n. iterations of chromaticty correction

%% For backward compatibilty
if (~isfield(DAoptions,'xmaxdas'))
    DAoptions.xmaxdas = 0.007;
end

if (~isfield(DAoptions,'xmindas'))
    DAoptions.xmindas = -0.015;
end

if (~isfield(DAoptions,'ymaxdas'))
    DAoptions.ymaxdas = 0.003;
end
   
if(~(isfield(DAoptions,'XmaxDA')))
    DAoptions.XmaxDA = 0.015;
end

if(~(isfield(DAoptions,'YmaxDA')))
    DAoptions.YmaxDA     = 0.004;
end


chrom_fams = LatticeOptData.chrom_fams; % list of sextupole families to be used for chromaticity correction

ACHRO       = LatticeOptData.ACHRO;
isdipole    = LatticeOptData.isdipole;

ACHRO  = setDVs(2, ACHRO,LatticeOptData, x);

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
       try
            [DA,~]=calcDA_fast(ACHRO,LatticeOptData.DAoptions,etax,betax,betay);
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
function LattStruct = cLatt(ACHRO,lattname,desc,LatticeOptData, MagnetStrengthLimits,varargin)
% Generates a structure for input to the lattice2LaFH function and
% optionally calculates and stores latice performance data.
%
%% Inputs
% Mandatory arguments
% ACHRO    : AT2 lattice cell array for one achromat. This may be a
%            structure following LatticeOptData.ACHRO (octupoles not included) 
%            or following LatticeOptData.ACHROGRD (octupoles included)
% lattname : latttice name
% desc     : lattice description
% LatticeOptData : structure with optimization data.
% MagnetStrengthLmits: structre with magnet strength Limits
%
% Optional flags
%
% 'basic'  : calculates basic lattice performnce data - atsummary
% 'DA'     : calculates dynamic apertures
% 'DAdist' : calculates DA distribution with errors
% 
%% Outputs
% LattStruc is a structure wih fields
% Lattice_Name = string wiht lattice name following naming convention
% Description : string with verbose description of how latice was developed
%
% ACHROMAT : AT2 lattice cell array for one achromat (an echo of the input)
%
% LattData.RINGGRD  : Full ring latice 6d enabled with same magnet strengths as
%            ACHROMAT. Octupoles are the same as in ACHRO (if they exist
%            there). Otehrwise they are zero.
% 
% LattData.lattMode : Lattice mode = 'a1,b1,b2,b3,c1,d1' (from LatticeOptData)
% LattData.Trb      : reverse bend angle [mrad]
% LattData.All_fams : cell array of string with the names of all families
% LattData.All_famsO := cell array of string with the names of all families including
%             octupoles
% LattData.XAll : strengths of all families (not including octupoles)
% LattData.XAllO : strengths of all families (not including octupoles)
% LattData.eqfam : cell array with the names of families listed in
%                  MagnetStrengthLimit corresponding to the All_fasmO for checking
%                  "challenge levels" with the chalevel.m function
% LattData.eqscal: scaling factors - useful to scald the dipole gradients before
%         comparison, when the dipole nedning angle has been changed to
%         compensate for the introduction of reverse bends.
%
% LattPerf.atsummary : output of atsummary run on RINGGRD 
%
%% Usage examples
% m4U_240316_b03_01_03_01=cLatt(ACHRO,'m4U_240316_b03_01_03_01',...
%        'MOGA_20240313T091640, Ind 12, Disp Matched, RB=3.49mrad',...
%         LatticeOptData,MagnetStrengthLimits);
% m4U_240331_b03_01_04_01=cLatt(rp.outputs.ACHROGRD,'m4U_240331_b03_01_04_01',...
%        'MOGA_20240330T062454, Ind 17, Disp Matched, Tunes Matched [56.20 19.28], RB=3.49mrad',...
%         LatticeOptData,MagnetStrengthLimits);
%

%% History
% PFT 2024/05 first version
% PFT 2024/06/03 : added info for magnet challenge level calculation
% PFT 2024/06/04 : updted output structure with new sub-fields
%
%% Input argument parsing

basicf  = any(strcmpi(varargin,'basic'));
DAf     = any(strcmpi(varargin,'DA'));
DAdistf = any(strcmpi(varargin,'DAdist'));

%% 
LattStruct=[];
LattStruct.Lattice_Name = lattname;
LattStruct.Description = desc;
LattStruct.ACHROMAT = ACHRO;
LattStruct.LattData.lattMode = LatticeOptData.lattMode;
LattStruct.LattData.Trb  = LatticeOptData.Trb;
LattStruct.LattData.All_fams = LatticeOptData.All_fams;
LattStruct.LattData.All_famsO = LatticeOptData.All_famsO;
LattStruct.LattData.eqfam = LatticeOptData.eqfam;
LattStruct.LattData.eqsca = LatticeOptData.eqsca;

if (length(ACHRO)==length(LatticeOptData.ACHRO))
    LattStruct.LattData.XAll  = getAllfams(2,ACHRO,LatticeOptData);
    LattStruct.LattData.XAllO = cat(2,LattStruct.LattData.XAll,[0.0 0.0 0.0]);
end
if (length(ACHRO)==length(LatticeOptData.ACHROGRD))
     LattStruct.LattData.XAll  = getAllfams(7,ACHRO,LatticeOptData);
     LattStruct.LattData.XAllO = getAllfamsO(7,ACHRO,LatticeOptData);
end
if ( (length(ACHRO)~=length(LatticeOptData.ACHRO))&&...
        (length(ACHRO)~=length(LatticeOptData.ACHROGRD)) )
        fprintf('%s Invalid input latice. Aborting...\n',datetime);
        return
end

RINGGRD = setAllfamsO(6,LatticeOptData.RINGGRD,LatticeOptData,LattStruct.LattData.XAllO);
LattStruct.LattData.RINGGRD=RINGGRD;
LattStruct.LattData.CLv=chalevel(LattStruct.LattData.XAllO,LatticeOptData.eqfam, ...
                         LatticeOptData.eqsca, MagnetStrengthLimits);


TolChrom     = LatticeOptData.DAoptions.TolChrom;% Chromaticity tolerances
Nitchro      = LatticeOptData.DAoptions.Nitchro; % Max n. iterations of chromaticity correction
chrom_fams   = LatticeOptData.chrom_fams;

[ACHR_zc, ~, ~]=fitchroit(ACHRO, chrom_fams, [0 0], Nitchro, TolChrom); 
LattStruct.LattData.ACHROMAT_ZC = ACHR_zc;

if (basicf)
    LattStruct.LattPerf.atsummary = atsummary(LattStruct.LattData.RINGGRD);
end

if (DAf)
    DAS_0 = calcDA(RINGGRD,LatticeOptData.DAoptions, 'dp', 0.00);
    DAS_p3 = calcDA(RINGGRD,LatticeOptData.DAoptions,'dp',+0.03);
    DAS_m3 = calcDA(RINGGRD,LatticeOptData.DAoptions,'dp',-0.03);
    LattStruct.LattPerf.DAS.DA_0=DAS_0;
    LattStruct.LattPerf.DAS.DA_p3=DAS_p3;   
    LattStruct.LattPerf.DAS.DA_m3=DAS_m3;
end




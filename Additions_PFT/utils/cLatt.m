function LattStruct = cLatt(ACHRO,lattname,desc,LatticeOptData)
% Generates a structure for input to the lattice2LaFH function
%% Inputs
% Mandatory argueents
% ACHRO    : AT2 lattice cell array for one achromat. This may be a
%            structure following LatticeOptData.ACHRO (octupoles not included) 
%            or following LatticeOptData.ACHROGRD (octupoles included)
% lattname : latttice name
% desc     : lattice description
% LatticeOptData : structure with optimization data.
%
%% Outputs
% LattStruc is a structure wih fields
% Lattice_Name = string wiht lattice name following naming convention
% Description : string with verbose desripton of how latice was developed
% ACHROMAT : AT2 lattice cell array for one achromat (an echo of the input)
% RINGGRD  : Full ring latice 6d enabled with same magnet strengths as
%            ACHROMAT. Octupoles are the same as in ACHRO (if they exist
%            there). Otehrwise they are zero.
% atsummary : output of atsummary run on RINGGRD 
% lattMode : Lattice mode = 'a1,b1,b2,b3,c1' (from LatticeOptData)
% Trb       = reverse bend angle [mrad]
% All_fams  = cell array of string wiht names of all families
% All_famsO = cell array of string wiht names of all families including
%             octupoles
% XAll = strengths of all families (not including octupoles)
% XAllO = strengths of all families (not including octupoles)
%
%% Usage examples
% m4U_240316_b03_01_03_01=cLatt(ACHRO,'m4U_240316_b03_01_03_01',...
%        'MOGA_20240313T091640, Ind 12, Disp Matched, RB=3.49mrad',...
%         LatticeOptData);
% m4U_240331_b03_01_04_01=cLatt(rp.outputs.ACHROGRD,'m4U_240331_b03_01_04_01',...
%        'MOGA_20240330T062454, Ind 17, Disp Matched, Tunes Matched [56.20 19.28], RB=3.49mrad',...
%         LatticeOptData);
%
LattStruct=[];
LattStruct.Lattice_Name = lattname;
LattStruct.Description = desc;
LattStruct.ACHROMAT = ACHRO;
LattStruct.lattMode = LatticeOptData.lattMode;
LattStruct.Trb  = LatticeOptData.Trb;
LattStruct.All_fams = LatticeOptData.All_fams;
LattStruct.All_famsO = LatticeOptData.All_famsO;

if (length(ACHRO)==length(LatticeOptData.ACHRO))
    LattStruct.XAll  = getAllfams(2,ACHRO,LatticeOptData);
    LattStruct.XAllO = cat(2,LattStruct.XAll,[0.0 0.0 0.0]);
else
    if (length(ACHRO)==length(LatticeOptData.ACHROGRD))
        LattStruct.XAll  = getAllfams(7,ACHRO,LatticeOptData);
        LattStruct.XAllO = getAllfamsO(7,ACHRO,LatticeOptData);
    else
        fprintf('%s Invalid input latice. Aborting...\n',datetime);
        return
    end
end
LattStruct.RINGGRD = setAllfamsO(6,LatticeOptData.RINGGRD,LatticeOptData,LattStruct.XAllO);
LattStruct.atsummary = atsummary(LattStruct.RINGGRD);
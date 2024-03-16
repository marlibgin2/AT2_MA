function LattStruct = cLatt(ACHRO,lattname,desc,LatticeOptData)
% Generates a structure for input to the lattice2LaFH function
%% input
% ACHRO    : AT2 lattice cell array for one achromat
% lattname : latttice name
% desc     : lattice descdription
% LatticeOptData : structure with optimization data.
%
%% Output
% LattStruc is a structure wih fields
% Lattice_Name = string wiht latice name following naming convention
% Description : string with verbose desripton of how latice was developed;
% ACHROMAT : AT2 lattice cell array for one achromat
% RINGGRD  : Full ring latice 6d enabled with same magnet strengths as
%            ACHROMAT
% atsummary : output of atsummary run on RINGGRD 
% lattMode : Lattice mode = 'a1,b1,b2,b3,c1' (from LatticeOptData)
% Trb  = reverse bend angle [mrad]
% All_fams = cell array of string wiht names of all families
% XAll = strengths of all families
%
%% Usage example
% m4U_240316_b03_01_03_01=cLatt(ACHRO,'m4U_20240316_b03_01_03_01',...
%        'MOGA_20240313T091640, Ind 12, Disp Matched', RB=3.49mrad,...
%         LatticeOptData);

LattStruct=[];
LattStruct.Lattice_Name = lattname;
LattStruct.Description = desc;
LattStruct.ACHROMAT = ACHRO;
XAll = getAllfams(2,ACHRO,LatticeOptData);
LattStruct.RINGGRD = setAllfams(6,LatticeOptData.RINGGRD,LatticeOptData,XAll);
LattStruct.atsummary = atsummary(LattStruct.RINGGRD);
LattStruct.lattMode = LatticeOptData.lattMode;
LattStruct.Trb  = LatticeOptData.Trb;
LattStruct.All_fams = LatticeOptData.All_fams;
LattStruct.XAll = getAllfams(2,ACHRO,LatticeOptData);

end
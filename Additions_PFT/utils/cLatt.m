function LattStruct = cLatt(ACHRO,lattname,desc,LatticeOptData)
% Generates a structure for input to the lattice2LaFH function
% input
% ACHRO    : AT2 lattice cell array
% lattname : latttice name
% desc     : lattice descdription
% LatticeOptData : structure with optimization data.
%
LattStruct=[];
LattStruct.Description = desc;
LattStruct.ACHROMAT = ACHRO;
LattStruct.RINGGRD = LatticeOptData.RINGGRD;
isdipole=LatticeOptData.isdipole;
LattStruct.ringpars = atsummary_fast(ACHRO,isdipole);
LattStruct.lattMode = LatticeOptData.lattMode;
LattStruct.Trb  = LatticeOptData.Trb;
LattStruct.XAll = getAllfams(2,ACHRO,LatticeOptData);
LattStruct.All_fams = LatticeOptData.All_fams;
LattStruct.RINGGRD = setAllfams(6,RINGGRD,LatticeOptData,XAll);
LattStruct.Lattice_Name = lattname;
end
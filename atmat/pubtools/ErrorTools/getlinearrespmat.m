function [ORM] = getlinearrespmat(RINGe)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

default_kick = 10e-6;

RINGData.Lattice = RINGe; 
RINGData.CavityFrequency = atGetRingProperties(RINGe,'rf_frequency'); 
RINGData.CavityHarmNumber = atGetRingProperties(RINGe,'harmonic_number');
CMData.HCMIndex = findcells(RINGe,'FamName','corrh'); 
CMData.VCMIndex = findcells(RINGe,'FamName','corrv');
CMData.HCMKicks = default_kick; 
CMData.VCMKicks = default_kick;
CMData.HCMCoupling = zeros(size(CMData.HCMIndex));
CMData.VCMCoupling = zeros(size(CMData.VCMIndex));
BPMData.BPMIndex = findcells(RINGe,'FamName','BPM');

RM = locoresponsematrix(RINGData,BPMData,CMData,'linear')./default_kick;

RM = mat2cell(RM, ...
    [ numel(BPMData.BPMIndex), numel(BPMData.BPMIndex)], ...
    [ numel(CMData.HCMIndex), numel(CMData.VCMIndex)]);

ORM.OrbHCor{1} = RM{1,1};
ORM.OrbHCor{3} = RM{2,1};
ORM.OrbVCor{1} = RM{1,2};
ORM.OrbVCor{3} = RM{2,2};
end
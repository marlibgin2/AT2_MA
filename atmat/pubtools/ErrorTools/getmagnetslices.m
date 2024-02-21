function MagnetIndices=getmagnetslices(RING)
%GETMAGNETSLICES Gets lattice indices for all magnet groups
% MagnetIndices = getmagnetslices(RING)
%
% INPUT
% 1. RING       AT2 lattice, with linked magnet elements sharing the same
%               MagNum field value.
%
% OUTPUT
% 1. MagnetIndices      {cell array of scalars} Each cell element contains
%                       the AT2 lattice indices for a specific magnet group.
% 
% NOTES
% 1. Rewritten version og getMagGroupsFromMagNum, to resolve indexing
% issues (last MagNum was skipped) 
%
% See also getMagGroupsFromMagNum


magnumind=findcells(RING,'MagNum');
magnumval=getcellstruct(RING,'MagNum',magnumind,1,1);

[~, ~, J] = unique(magnumval);
MagnetIndices = cell(numel(unique(magnumval)),1);

for n = 1:numel(J)
    MagnetIndices{J(n)} = cat(2,MagnetIndices{J(n)},magnumind(n)); 
end
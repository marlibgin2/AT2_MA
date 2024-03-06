function RING = markSlicedMagnets(RING,varargin)
% MARKSLICEDMAGNETS adds a MagNum field to all sliced elements in a lattice
% 
%  NEWRING = markSlicedMagnets(RING,StartingID)
%
%  INPUT
%  1. RING          AT2 lattice
%  2. StartingID    {scalar, optional} Starting index to use. Added ID
%                   numbers will then increase by 1 for each magnet.
%
%  OUTPUT
%  2. NEWRING       AT2 lattice, with added MagNum field for sliced magnets
%
%  NOTES
%  1. The function assumes that any elements with magnetic fields (e.g.
%  PolynomA, PolynomB or K) with the same FamName that are in direct
%  contact, i.e. with no non-zero length elements in between them,
%  constitute a sliced magnet. All such slices are assigned their own
%  unique MagNum ID for further use by error functions.
%  2. Remember that any subsequent duplication of a lattice, such as via
%  achromat2ring, will not update the MagNum field for duplicated elements.
%  
%  See also achromat2ring, getMagGroupsFromMagNum

if nargin > 1
    MagNum = varargin{1};
else
    MagNum = 1;
end

% Helper functions
    function y = hasField(x)
        y = (isfield(RING{x},'PolynomB') && any(RING{x}.PolynomB ~= 0)) || ...
            (isfield(RING{x},'PolynomA') && any(RING{x}.PolynomA ~= 0)) || ...
            (isfield(RING{x},'K') && RING{x}.K ~= 0);
    end

    function y = previousThickElement(x)
        for n = (x-1):-1:(x - 1 -numel(RING))
            if RING{mod(n-1,numel(RING))+1}.Length ~= 0
                y = n;
                return;
            end
        end
    end
    

% Main loop
for k = 1:numel(RING)


    % If this element has a magnetic field, check the preceding thick
    % element if that also has a field. As the assumption is neighbouring
    % elements with fields are sliced magnets, use the preceding MagNum in
    % that case. Or assign a new one, if there was none.
    if hasField(k)
        p = previousThickElement(k);
        if strcmpi(RING{k}.FamName,RING{p}.FamName)
            RING{k}.MagNum = RING{p}.MagNum;
        else
            RING{k}.MagNum = MagNum;
            MagNum = MagNum + 1;
        end
    end
end
end

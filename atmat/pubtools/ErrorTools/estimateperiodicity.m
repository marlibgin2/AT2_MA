function varargout = estimateperiodicity(RING)
% ESTIMATEPERIODICITY determines periodicity based on the horizontal beta
%  [N, index] = estimateperiodicity(RING)
%
% INPUT
% 1. RING       {cell array of structs} AT2 lattice variable
%
% OUTPUT
% 1. N          {scalar, optional} Estimated ring periodicity
% 2. index      {cell array} lattice element indices, sorted into a cell
%               array of length N depending on which achromat they belong to
%
% NOTES AND LIMITATIONS
% 1. Tested for MAX IV lattices (R1 and R3)
% 2. Relies on FFT of the hor. beta function and computes the periodicity
%    based on the lowest occuring frequency. Therefore, lattices with
%    superperiods will require code adaptation (e.g. Diamond)!
% 3. Sorting of elements into achromats is based on s-coordinate. Therefore
%    thin elements at achromat boundaries may be allocated to the wrong 
%    achromat.



LD = linopt(RING,1e-4,1:numel(RING)+1);

S = cat(2,LD.SPos);
[~, ai] = unique(S);
S = S(ai);
BETA = cat(1,LD.beta); BETA = BETA(ai,1);

% Interpolate to get even sampling; unclear on whether this is needed
nbrPoints = (1*2*3*5)^3;    % Interpolation points. Might be overkill.
ds = S(end)/(nbrPoints-1);
Si = 0:ds:S(end);
BETAi = interp1(S,BETA,Si);

% Do an FFT, grab the first peak on either side and calculate the
% periodicity based on this.
FBETA = fft(BETAi - mean(BETAi));
[ amps, pks ] = findpeaks(abs(FBETA));
[~, k] = sort(amps(1:end));
k = flip(k);
N = (nbrPoints - max(pks(k(1:2))) + min(pks(k(1:2))))/2;


if nargout > 0
    varargout{1} = N;
end

if nargout > 1
    % Get the achromat element indices
    S = cat(2,LD.SPos);
    C = S(end);
    achrL = C/N;

    achrIndices = cell(N,1);
    for n = 1:N
        startPos = (n-1)*achrL;
        endPos = n*achrL;
        Index = S >= startPos & S < endPos;
        achrIndices{n} = find(Index);
    end
    varargout{2} = achrIndices;
end

end
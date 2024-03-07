function Xout = trunc_randn(varargin)
%TRUNC_RANDN generate a normal gaussian distribution truncated at N sigma
% Xout = trunc_randn(Nel0, Nsigma)
%
% INPUT
% 1. Nel0   - Length of output vector
% 2. Nsigma - Number of sigma above which to truncate {default: 2}
%
% OUTPUT
% 1. Xout   - Vector containing the output
%
% NOTES
% 1. Used in manufacturing scenarios where one can assume random errors and
% a Q/A process able to catch and correct deviations larger than a certain
% amount.
%
% See also randn


Nel0 = varargin{1};

if nargin<2
    Nsigma=2;
else
    Nsigma=varargin{2};
end
Nel = Nel0;
% Xout = [];
% while Nel>0
%     X = randn(Nel,1);
%     X = X(abs(X)<Nsigma);
%     Xout = [Xout; X]; %#ok<AGROW> 
%     Nel   = Nel - max(size(X));
% end

% Rewritten to use recursive calls. The above appears to have some sort of
% bug that occasionally drops an element
X = randn(Nel,1);
outOfBounds = abs(X) >= Nsigma;
if any(outOfBounds)
    X(outOfBounds) = trunc_randn(sum(outOfBounds),Nsigma);
end
Xout = X;

end


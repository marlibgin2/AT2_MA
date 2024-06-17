function varargout = mcf(RING,dp0)
%MCF momentum compaction factor
% MCF(RING) calculates momentum compaction factor of RING
%
% MCF(RING,DPP) computes the momentum compaction for off-momentum DPP
%
% [A1, A2, ...] = MCF(RING, ...) computes higher order mom. compaction

if nargin < 2, dp0=0; end
mcforder = nargout;

% Establish the momentum shift vector. Numerically the range has to be
% wider the more orders are requested, otherwise the higher coefficients
% will be fitted to noise. This is dealt with via a somewhat empirical
% recipe.
ddp = 1e-6 * sqrt(10)^mcforder; 
ddp = (-0.5:1/mcforder:0.5)*ddp;

% Build initial particle positions that begin on the closed orbit
% transversally, but with energy shifts and with cT = 0
X0 = nan(6,mcforder);
for n = 1:(mcforder+1)
    X0(1:4,n) = findorbit4(RING,dp0 + ddp(n));
    X0(5,n) = dp0 + ddp(n);
    X0(6,n) = 0;
end

% Some orbits may not have been found. Assign nan to the outputs and
% return. Leave it as an exercise for the user whether he/she wishes to
% decrease the order or offset.
if any(isnan(X0(:)))
    varargout = num2cell(nan(nargout,1));
    return;
end

% Track X0 for 1 turn for all energy shifts to determine the path
% difference (6th canonical variable)
T = ringpass(RING,X0);
dL = T(6,:)';

% Get circumference
L = findspos(RING,length(RING)+1);

% Do a polynomial fit (rescales internally for robustness)...
[phat,~,mu] = polyfit(ddp',dL,mcforder);
p = rescale(phat,mu); ..., then transform back

% Divide the obtained coefficients with the circumference as usual to get
% the momentum compaction coefficients
a = fliplr(p)./L;

% Finally, ditch the constant (should be zero) and return the answer
varargout = num2cell(a(2:end));
end


%% HELPER FUNCTION
function p = rescale(phat,mu)

phat = flip(phat); % Flip order to [p0, ..., pn]
n = numel(phat)-1;
p = zeros(size(phat));
for i = 0:n
    for k = i:n
        p(i+1) = p(i+1) + nchoosek(k, k-i) * phat(k+1)/mu(2)^k * (-mu(1))^(k-i);
    end
end
p = flip(p); % Back to original order [pn, ..., p0]
end


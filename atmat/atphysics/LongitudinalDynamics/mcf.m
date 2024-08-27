function varargout = mcf(RING,varargin)
%MCF momentum compaction factor
% MCF(RING) calculates momentum compaction factor of RING
%
% MCF(RING,DPP) computes the momentum compaction for off-momentum DPP
% MCF(RING,...,'Display') plots dL(dpp) along with fitted mcf polynomial.
%
% [A1, A2, ...] = MCF(RING, ...) computes higher order mom. compaction

% if nargin < 2, dp0=0; end
mcforder = nargout;

dp_range = 0.1;         % Max p-p energy excursion to consider.
dp0 = 0;                % Center energy around which to fit mom. cmp.
max_dp_step = 0.01;     % Limit on how sparse an energy sampling is allowed
DisplayFlag = false;    % Whether to show diagnostic information

for n = nargin-1:-1:1
    if ischar(varargin{n})
        switch lower(varargin{n})
            case 'display'
                DisplayFlag = 1;
        end
    elseif isnumeric(varargin{n})
        dp0 = varargin{n};
    end
end

% Establish the momentum shift vector. Numerically the range has to be
% wider the more orders are requested, otherwise the higher coefficients
% will be fitted to noise. This is dealt with via a somewhat empirical
% recipe.
if mcforder == 1
    dp_range = 1e-6;     % This is used to match the behaviour of the baseline mcf function, i.e. preserve previous behaviour
end
ddp = (-0.5:min(0.5/mcforder,max_dp_step/dp_range):0.5)*dp_range;

% Build initial particle positions that begin on the closed orbit
% transversally, but with energy shifts and with cT = 0
X0 = nan(6,numel(ddp));
for n = 1:numel(ddp)
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

% Filter away the non-ok orbits and energies
iok = ~sum(isnan(X0));
X0 = X0(:,iok);
ddp = ddp(iok);

% Track X0 for 1 turn for all energy shifts to determine the path
% difference (6th canonical variable)
T = ringpass(atdisable_6d(RING),X0);
dL = T(6,:)';

% Get circumference
L = findspos(RING,length(RING)+1);

% Do a polynomial fit (rescales internally for robustness)...
[phat,~,mu] = polyfit(ddp',dL,mcforder);
p = rescale(phat,mu); ..., then transform back

% If requested, plot dL(ddp) along with the fit
if DisplayFlag
    F = [];
    for n = 0:mcforder
        F = cat(1,ddp.^n,F);
    end
    
    figure; 
    
    subplot(2,1,1);
    plot(ddp,dL,'xb','DisplayName','cT from ringpass'); hold on; plot(ddp,p*F,'-k','DisplayName','Polynomial fit');
    xlabel('\delta');
    ylabel('cT [m]');
    legend('Location','northwest');
    grid on;

    subplot(2,1,2);
    plot(ddp,dL' - p*F);
    xlabel('\delta');
    ylabel('Residual [m]');
end


% Divide the obtained coefficients with the circumference as usual to get
% the momentum compaction coefficients
a = fliplr(p)./L;

% Finally, ditch the constant (should be zero) and return the answer
varargout = num2cell(a(2:nargout+1));
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

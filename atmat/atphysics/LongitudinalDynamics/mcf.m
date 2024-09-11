function varargout = mcf(RING,varargin)
%MCF momentum compaction factor
% MCF(RING) calculates momentum compaction factor of RING
%
% MCF(RING,DPP) computes the momentum compaction for off-momentum DPP
%
% MCF(RING,...,'Display') plots dL(dpp) along with fitted mcf polynomial.
%
% [A1, A2, ...] = MCF(RING, ...) computes higher order mom. compaction
%
% MCF(RING,...,'Range',dp_range) specifies the range over which the higher
%                                order mom. comp. factors should be
%                                calculated. Default is 10% (+-5% around
%                                on-momentum particle) 
%
% NOTES
% 1. If only linear momentum compaction is requested dp_range is set to
%    1e-6. This is for reasons of backwards compatibility. In addition 
% 2. For lattices with a limited momentum acceptance it may be necessary to
%    use the 'Range' argument to reduce the momentum span when calculating
%    higher order mcf. Use momentum_aperture_at to determine suitable
%    boundaries. 
%
% IMPORTANT!!!
% MCF gives a wrong result with 6-d rings. The RING should be set to 4d.
%
% See also: ATDISABLE_6D, CHECK_6D, MOMENTUM_APERTURE_AT

mcforder = max(nargout,1);  % Ensure there is always at least one output, to ensure backwards compatibility

%% DEFAULT PARAMETERS

% Max p-p energy excursion to consider. In order to preserve previous
% behaviour, unless otherwise requested further down, if only the linear
% mcf is requested the energy excursions is limited to +-1 ppm.
if mcforder == 1
    dp_range = 1e-6;     
else
    dp_range = 0.1;     
end

dp0 = 0;                % Center energy around which to fit mom. cmp.
max_dp_step = 0.01;     % Limit on how sparse an energy sampling is allowed
DisplayFlag = false;    % Whether to show diagnostic information

%% INPUT HANDLING
for n = nargin-1:-1:1
    if ischar(varargin{n})
        switch lower(varargin{n})
            case 'display'
                DisplayFlag = 1;
            case 'range'
                dp_range = varargin{n+1};
        end
    elseif isnumeric(varargin{n})
        if n == 1
            dp0 = varargin{n};
        end
    end
end

%% COMPUTATIONS
% Establish the momentum shift vector. Numerically the range has to be
% wider the more orders are requested, otherwise the higher coefficients
% will be fitted to noise. This is dealt with via a somewhat empirical
% recipe.
ddp = (-0.5:min(0.5/(mcforder),max_dp_step/dp_range):0.5)*dp_range;

% Build initial particle positions that begin on the closed orbit
% transversally, but with energy shifts and with cT = 0
X0 = nan(6,numel(ddp));
for n = 1:numel(ddp)
    X0(1:4,n) = findorbit4(RING,dp0 + ddp(n),'strict', -1);
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
T = ringpass(RING,X0,'KeepLattice');
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
varargout = num2cell(a(2:mcforder+1));
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

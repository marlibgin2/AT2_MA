function tunes=calcTune(RING,Rin,varargin)
% Calculates betatron tunes for a set of initial coordinates
% based on atnuampl:
% 
% Differences wrt atnuampl
%   Tracking is done in absolute coordinates, not wrt the closed orbit.
%   added ability ability to deal with a grid of points 
%   added ability to do off-energy tracking.
%
%% Inputs
% Mandatory arguments
%
% RING: AT2.0 lattice structure
% Rin: (6XN) matrix of initial coordinates 
%
% Optional arguments
%   nturns :# of turns for tracking for each set, default = 256.
%   nsets  :# of sets of consecutive turns for which tunes are to be
%           determined, default = 1
%   minampx: minimum initial horizontal particle amplitude, default = 30 µm.
%   minampy: minimum initial vertical particle amplitude, default = 30 µm.
%   method: specify the method for tune determination
%        1: Highest peak in fft
%        2: Interpolation on fft results
%        3: Windowing + interpolation (default)
%        4: NAFF
%% Optional flags
% verbose : turns on verbose output
%
%% Outputs
% tunes structure with fields
% tunes.inputs echoes the input parameters
%   tunes.inputs.RING
%   tunes.inputs.Rin
%   tunes.inputs.nturns
%   tunes.inputs.nsets
%   tunes.inputs.minampx
%   tunes.inputs.minampy
%   tunes.inputs.method
%
% tunes.outputs.Qx  : (Nxnsets) array of Horizontal tunes 
% tunes.outputs.Qy  : (Nxnsets) array vertical tunes
% tunes.outputs.lost: (Nx1) array of logicals (1= particle was lost)
% tunes.outputs.Rout: Particle coordinates at the end of tracking. 
%
% NOTES
% 1. Be aware that when looking at tuneshifts due to an oscillation
%    amplitude in the orthonormal plane there is a risk of the ADTS
%    calculation picking the wrong peak once the amplitudes become large
%    enough due to betatron coupling. This is dealt with in the NAFF method
%    by imposing a threshold on tune amplitude shifts from one amplitude to
%    the next, before picking the largest amplitude peak.
%
% See also findtune, calcnaff
%% Usage examples
% tunes = calcTune(RING,[0.001 0 0 0 0 0 0]','nturns', 128, 'method', 3);
% tunes = calcTune(RING,Rin, 'nturns',1024, 'nsets',2, 'method', 4);


%% History
% PFT 2024/04/28 First version, based on atnuampl with modifications by
%                M.Sjöström
% PFT 2024/05/03 Added optional flag for enabling/suppressing verbose output
% PFT 2024/05/04 Added tune calcualtion for two consectuvie sets of turns 
%                - Useful for tune  diffusion map  calculations

%% Input argument parsing
nturns   = getoption(varargin,'nturns',128);
nsets    = getoption(varargin,'nsets', 1);
minampx  = getoption(varargin,'minampx',30E-6);
minampy  = getoption(varargin,'minampy',30E-6);
method   = getoption(varargin,'method',3);
verbosef = any(strcmpi(varargin,'verbose'));

tunes.inputs.RING    = RING;
tunes.inputs.nturns  = nturns;
tunes.inputs.minampx = minampx;
tunes.inputs.minampy = minampy;
tunes.inputs.method  = method;

%% Preamble
npt = size(Rin,2);

if method == 4    % Turn adjustment, recommended for NAFF
    nturns = 2^(log2(nturns));
    nturns = nturns + 6-mod(nturns,6);
end

[~,nbper]=atenergy(RING);
[lindata,fractune0]=atlinopt(RING,0.0,1:length(RING)+1);
tune0=nbper*lindata(end).mu/2/pi;
offs=repmat([nbper -nbper],1,nsets);

Rin(1,(Rin(1,:)==0))=minampx;
Rin(3,(Rin(3,:)==0))=minampy;

%% Tracks and calculate tunes

[p1, lost]=ringpass(RING,Rin,nturns*nsets);
Rout = p1(:,1:nturns:end);

tunetrack  = nan(npt,2*nsets);
tunetrackm = [];
for j=1:nsets
  if (method==4)  
    for n=1:npt
        if (lost(n))
        % Particle was lost during tracking, don't run NAFF, assign NaN
        % and continue 
            tunetrack(n,1:2*nsets) = nan(1,2*nsets);
            continue;
        end
        particleTurns = (n+npt*nturns*(j-1)):npt:(npt*nturns+npt*nturns*(j-1)); 
    
        [nux, amplx, ~] = calcnaff(p1(1,particleTurns),p1(2,particleTurns));
        [nuy, amply, ~] = calcnaff(p1(3,particleTurns),p1(4,particleTurns));
        
        % If NAFF fails for any reason assign NaN and continue
        if any(isnan([nux;nuy]))
            tunetrack(n,1+nsets*(j-1):2+nsets*(j-1)) = nan(1,2);
            continue;
        end

        if n > 1
        % If there is a tune calculated at a lower amplitude, first
        % filter the tune peaks from NAFF to ignore frequencies
        % differing too much from the previous amplitude. The threshold
        % has been arbitrarily set to 0.1.
        % This is useful primarily when looking at the horizontal tune
        % shift from a vertical kick, or vice versa, in a lattice with
        % a large amount of coupling.
          dnux =  abs(tunetrack(n-1,1+2*(j-1)) - abs(nux)./(2*pi));
          dnuy =  abs(tunetrack(n-1,2+2*(j-1)) - abs(nuy)./(2*pi));
          ix = dnux < 0.1; nux = nux(ix); amplx = amplx(ix);
          iy = dnuy < 0.1; nuy = nuy(iy); amply = amply(iy);
       end
       [~, i] = max(amplx);        % Identify the dominant frequency peak
       nux = abs(nux(i))/(2*pi);   % Re-normalize to get the tune. Note the sign is ignored.

       [~, i] = max(amply);
       nuy = abs(nuy(i))/(2*pi);
        
       tunetrack(n,1+nsets*(j-1):2+nsets*(j-1)) = [nux,nuy];
    end
  else

    if (verbosef)
        tunetrackj=[findtune(reshape(p1(1,1+npt*nturns*(j-1):npt*nturns*j),npt,nturns)','verbose',method);...
                   findtune(reshape(p1(3,1+npt*nturns*(j-1):npt*nturns*j),npt,nturns)','verbose',method)]';
    else
        tunetrackj=[findtune(reshape(p1(1,1+npt*nturns*(j-1):npt*nturns*j),npt,nturns)',method);...
                   findtune(reshape(p1(3,1+npt*nturns*(j-1):npt*nturns*j),npt,nturns)',method)]';
    end
    tunetrackm = [tunetrackm tunetrackj];
  end
end

if (not(method==4))
    tunetrack=tunetrackm;
end

[~,k]=min([repmat(fractune0,1,nsets)-tunetrack(1,:); 1-repmat(fractune0,1,nsets)-tunetrack(1,:)]);
np=offs(k);
offset=round(repmat(tune0,1,nsets)-np.*tunetrack(1,:));
tunetrack=np(ones(npt,1),:).*tunetrack + offset(ones(npt,1),:);

%% Collects output structure data
tunes.outputs.Qx   = tunetrack(:,1:nsets:2*nsets-1);
tunes.outputs.Qy   = tunetrack(:,2:nsets:2*nsets);
tunes.outputs.lost = lost';
tunes.outputs.Rout = Rout;



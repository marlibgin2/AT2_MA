function varargout=atnuampl(ring,ampl,xy,varargin)
%ATNUAMPL	computes tune shift with amplitude
%[NUX,NUY]=ATNUAMPL(RING,AMPLITUDE)
%[NUX,NUY]=ATNUAMPL(RING,AMPLITUDE,1)
%
%	Computes tunes for the specified horizontal amplitudes
%
%[NUX,NUY]=ATNUAMPL(RING,AMPLITUDE,3)
%
%	Computes tunes for the specified vertical amplitudes
%
%ATNUAMPL(...)
%   Plots the computed tunes in the current axes
%
%ATNUAMPL(...,Name,Value)
%   Uses additional options specified by one or more Name,Value pairs.
%   Possible values are:
%       orbit:  initial closed orbit
%       nturns: specify the number of turns for tracking (default 256)
%       minamp: specify the min. initial particle amplitude (default 30 µm)
%       method: specify the method for tune determination
%               1: Highest peak in fft
%               2: Interpolation on fft results
%               3: Windowing + interpolation (default)
%               4: NAFF
%   Other options are transmitted to the plot function
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


lab={'x^2','p_x^2','z^2','p_z^2'};
if nargin < 3, xy=1; end
[nturns,varargs]=getoption(varargin,'nturns',256);
[method,varargs]=getoption(varargs,'method',3);
[orbit,varargs]=getoption(varargs,'orbit',[]);
[minAmp,varargs]=getoption(varargs,'minamp',30e-6);

if method == 4    % Turn adjustment, recommended for NAFF
    nturns = 2^(log2(nturns));
    nturns = nturns + 6-mod(nturns,6);
end

if ~isempty(varargs) && isnumeric(varargs{1})	% ATNUAMPL(RING,AMPLITUDE,XZ,ORBIT)
    orbit = varargs{1};
    varargs(1)=[];
end

if isempty(orbit)
    if check_radiation(ring)
        orbit=findorbit6(ring);
        dp=orbit(5);
    else
        dp=0.0;
        [~, orbit]=findorbit4(ring, dp);
    end
end

[~,nbper]=atenergy(ring);
[lindata,fractune0]=atlinopt(ring,dp,1:length(ring)+1, 'orbit', orbit);
tune0=nbper*lindata(end).mu/2/pi;
offs=[nbper -nbper];
siza=size(ampl);
nampl=prod(siza);
p0=repmat(minAmp*[1;0;1;0;0;0], 1,nampl); % specifies minimum amplitude
%p0(xy,:)=max(p0(xy,:),ampl(:)');
ampl(ampl==0)=minAmp;
p0(xy,:)=ampl;
p0=p0+orbit(:,ones(1,nampl));
p1=ringpass(ring,p0,nturns)-orbit(:,ones(1,nampl*nturns));
if method == 4
    tunetrack = nan(numel(ampl),2);
    for n = 1:numel(ampl)
        particleTurns = n:nampl:(numel(ampl)*nturns);

        if any(isnan(p1(1:4,particleTurns)),'all')
            % Particle was lost during tracking, don't run NAFF, assign NaN
            % and continue 
            tunetrack(n,1:2) = nan(1,2);
            continue;
        end

        [nux, amplx, ~] = calcnaff(p1(1,particleTurns),p1(2,particleTurns));
        [nuy, amply, ~] = calcnaff(p1(3,particleTurns),p1(4,particleTurns));
        
        % If NAFF fails for any reason assign NaN and continue
        if any(isnan([nux;nuy]))
            tunetrack(n,1:2) = nan(1,2);
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
            dnux =  abs(tunetrack(n-1,1) - abs(nux)./(2*pi));
            dnuy =  abs(tunetrack(n-1,2) - abs(nuy)./(2*pi));
            ix = dnux < 0.1; nux = nux(ix); amplx = amplx(ix);
            iy = dnuy < 0.1; nuy = nuy(iy); amply = amply(iy);
        end
        [~, i] = max(amplx);        % Identify the dominant frequency peak
        nux = abs(nux(i))/(2*pi);   % Re-normalize to get the tune. Note the sign is ignored.

        [~, i] = max(amply);
        nuy = abs(nuy(i))/(2*pi);
        
        tunetrack(n,1:2) = [nux,nuy];
    end
     
else
    tunetrack=[findtune(reshape(p1(1,:),nampl,nturns)',method,'verbose',false);...
        findtune(reshape(p1(3,:),nampl,nturns)',method,'verbose',false)]';
end
[~,k]=min([fractune0-tunetrack(1,:); 1-fractune0-tunetrack(1,:)]);
np=offs(k);
offset=round(tune0-np.*tunetrack(1,:));
tunetrack=np(ones(nampl,1),:).*tunetrack + offset(ones(nampl,1),:);
if nargout > 0
    varargout={reshape(tunetrack(:,1),siza),reshape(tunetrack(:,2),siza)};
else
    inttunes=floor(tune0);
    plot((ampl.*ampl)',tunetrack-inttunes(ones(nampl,1),:),'o-',varargs{:});
    legend('\nu_x','\nu_z');
    xlabel(lab{xy});
    ylabel('\nu');
    grid on
end


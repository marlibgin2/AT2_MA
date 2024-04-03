function varargout=atnuampl(ring,ampl,xz,varargin)
%ATNUAMPL	computes tune shift with amplitude
%[NUX,NUZ]=ATNUAMPL(RING,AMPLITUDE)
%[NUX,NUZ]=ATNUAMPL(RING,AMPLITUDE,1)
%
%	Computes tunes for the specified horizontal amplitudes
%
%[NUX,NUZ]=ATNUAMPL(RING,AMPLITUDE,3)
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
%       method: specify the method for tune determination
%               1: Highest peak in fft
%               2: Interpolation on fft results
%               3: Windowing + interpolation (default)
%               4: NAFF
%   Other options are transmitted to the plot function
%
% See also findtune, calcnaff


lab={'x^2','p_x^2','z^2','p_z^2'};
if nargin < 3, xz=1; end
[nturns,varargs]=getoption(varargin,'nturns',256);
[method,varargs]=getoption(varargs,'method',3);
[orbit,varargs]=getoption(varargs,'orbit',[]);

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
p0=repmat(0.00003*[1;0;1;0;0;0], 1,nampl); % 30 microns minimum amplitude
%p0(xz,:)=max(p0(xz,:),ampl(:)');
ampl(ampl==0)=0.00003;
p0(xz,:)=ampl;
p0=p0+orbit(:,ones(1,nampl));
p1=ringpass(ring,p0,nturns)-orbit(:,ones(1,nampl*nturns));
if method == 4
    tunetrack = nan(numel(ampl),2);
    for n = 1:numel(ampl)
        particleTurns = n:nampl:(numel(ampl)*nturns);
        [nux, amplx, ~] = calcnaff(p1(1,particleTurns),p1(2,particleTurns));
        [~, i] = max(amplx);        % Identify the dominant frequency peak
        nux = abs(nux(i))/(2*pi);   % Re-normalize to get the tune. Note the sign is ignored.

        [nuy, amply, ~] = calcnaff(p1(3,particleTurns),p1(4,particleTurns));
        [~, i] = max(amply);
        nuy = abs(nuy(i))/(2*pi);
    
        tunetrack(n,1:2) = [nux,nuy];
    end
     
else
    tunetrack=[findtune(reshape(p1(1,:),nampl,nturns)',method);...
        findtune(reshape(p1(3,:),nampl,nturns)',method)]';
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
    xlabel(lab{xz});
    ylabel('\nu');
    grid on
end

function [tune, zw] = naff2(z, guess, tol)
% numerical analysis of fundamental frequencies
% [tune, zw] = naff(z, top)
% z : complex(samples, bpms)
% tol (default = 1e-15) : numerical tolerance 

if nargin < 2
    guess = NaN;
end

if nargin < 3
    tol = 1e-15;
end
    
maxn = length(z);
mft = ceil(log(maxn)/log(2));
npoint = 2 .^ mft;

% hann window
st = 2*pi/maxn;
hannw = 1+cos(st*((1:maxn)-maxn/2));
z = z.*hannw;

% fft take highest bin as tune guess
zsing = fft(z, npoint)./npoint;

% only look inside guess window if provided
if ~isnan(guess)
    wy = floor(guess*length(zsing));
    wy = wy(1):wy(2);
    wy(wy < 1) = 1;
    wy(wy > length(zsing)) = length(zsing);
    [ftmax, nftmax] = max(abs(zsing(wy)));
    nftmax = nftmax + wy(1) - 1;
else
    [ftmax, nftmax] = max(abs(zsing));
end

tuneguess = (nftmax-1)./npoint;

deltat = 1 ./ npoint;
tuneguess(tuneguess > 0.5) = tuneguess(tuneguess > 0.5)-1;

% tune1 is left bound of tune
tune1 = tuneguess - deltat;

% tune at maximum of absolute value of sum of phasor * signal

iter = 0;

function m = phasmag(v)
    ztune = exp(-2.*pi.*1j.*v);
    zf = z(maxn);
    for n = maxn-1:-1:1
        zf = zf*ztune+z(n);
    end
    m = -abs(zf);
    iter = iter + 1;
end

% search without using derivatives

o = fminbnd('defaults');
o.TolX = tol;

[tune, zw, flag] = fminbnd(@phasmag, tune1, tune1 + deltat * 2, o);
    
% problems?

if flag < 0
    tune = tuneguess;
    zw = ftmax;
end

% disp(iter);
zw = -zw ./ npoint;

ztune = exp(-2.*pi.*1j.*tune);
zf = z(maxn);
for nn = maxn-1:-1:1
    zf = zf*ztune+z(nn);
end
zw = zf ./ npoint;

end
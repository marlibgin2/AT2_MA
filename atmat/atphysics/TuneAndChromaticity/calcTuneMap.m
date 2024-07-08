function tunemap=calcTuneMap(varargin)
% Calculates and plots betatron tunes maps as a function of position 
% or energy.
% Calculates an plots frequency maps and tune diffusion maps
%
% This is a high level wrapper to the function calcTune. 
% Tracking is 4d + fixed energy deviation
%   
%% Inputs 
% Mandatory arguments
% ACHRO: AT2 lattice array. 
% TMoptions: structure with fields (if TMoptions is empty, defaults are
%                                   used. any specific fields may be
%                                   changed by entering an optional
%                                   argument)
%  TMoptions.mode: (default is 'x'):
%  'x'      calculates only tunes vs horizontal position at minimum
%           vertical position (minampy)
%  'y'      calculates only tunes vs vertical position at minimum
%           horizontal position (minampx)
%  'xy'     calculates both x and y as described above
%  'gridxy' calculates tunes on a grid of points in (x,y) plane.
%  'gridxdp'calculates tunes on a grid of points in (x,dp) plane
%  'gridydp'calculates tunes on a grid of points in (x,dp) plane
%  'difxy'  calculates tune diffusion map in the xy plane 
%               tune calculation method is 4 (NAFF). 
%  'difxdp' calculates tune diffusion map in the xdp plane
%               tune calculation method is 4 (NAFF). 
%  'difydp' calculates tune diffusion map in the ydp plane 
%               tune calculation method is 4 (NAFF). 
%  'chrom' calculates chromatic tune footprint.
%  TMoptions.minampx: minimum absolute value of amplitude in horizontal direction,
%            default=30 microm
%  TMoptions.minampy: minimum absolute value of amplitude in horizontal direction,
%            default=30 microm
%  TMoptions.xmax: max horizontal amplitude [m], default = 0.005 
%  TMoptions.xmin: min horizontal amplitude [m], default = 0.0 
%       
%  TMoptions.ymax: max vertical amplitude [m], default = 0.004
%  TMoptions.ymin: min vertical amplitude [m], default = 0.0
% 
%  TMoptions.dpmax: max energy deviation for chromatic tune footprint,default = +0.03
%  TMoptions.dpmin: min energy deviation for chromatic tune footprint,default = -0.03
%
%  TMoptions.dp: initial energy deviation, default = 0.0
%
%  TMoptions.npx: number of points along horizontal direction; default = 11
%  TMoptions.npy: number of points along vertical direction: default = 11
%  TMoptions.npd: number of points along energy axis: default = 11
%  TMoptions.nturns : number of turns, default = 128
%
%  TMoptions.method 1: Highest peak in fft
%                   2: Interpolation on fft results
%                   3: Windowing + interpolation (default)
%                   4: NAFF (this is always the method in case the mode is "diff"
%
%  TMoptions.smooth  : if true, selects smooth mode grid calculations - note that this means 
%           computations are not parallelized. This is only relevant for
%           method=4 (NAFF), deafult = false
%
% Optional arguments
% Any of the fields in the TMoptions structure
% desc: descriptive string
%
% Optional ploting options 
% plotmode: 'abs' : plots full tune value (inc. integer part) (default)
%           'rel' : plots tune variations with respect to small amplitude
%                   tunes.
% plottype: 
%   'x'       : ADTS along horizontal axis (default)
%   'y'       : ADTS along vertical axis
%   'xy'      : both ADTS
%   'td'      : same as 'xy', but on a tune diagram
%   'gridtd'  : grid of points on a tune diagram
%   'gridxy'  : grid of points as a series of lines vs 
%               x for a chosen set of y's
%   'gridyx'  : grid of points as a series of lines vs 
%               y for a chosen set of x's
%   'gridxdp' : grid of points as a series of lines vs 
%               x for a chosen set of dp's
%   'gridydp' : grid of points as a series of lines vs 
%               y for a chosen set of dp's
%   'griddpx' : grid of points as a series of lines vs 
%               dp for a chosen set of x's
%   'griddpy' : grid of points as a series of lines vs 
%               dp for a chosen set of y's
%   'difxy' : tune diffusion plot on xy plane
%   'difxdp': tune diffusion plot on xdp plane
%   'difydp': tune diffusion plot on ydp plane
%   'fmxy'  : xy tune diffusion plot on tune diagram (frequency map)
%   'fmxdp' : xdp tune diffusion plot on tune diagram (frequency map) 
%   'fmydp' : xdp tune diffusion plot on tune diagram (frequency map)
%   'chro'  : tunes vs energy deviation
%   'chrotd': chromatic tune footprint on a tune diagram
%
%
% resorder: resonance order for tune diagram, default = 5
% qxrange=[qxmin qymin]: horizontal plot range in tune diagram,  default =[0 1]
% qyrange=[qymin qymax]: vertical plot range in tune diagram, default= [0 1]
% xminplot: minimum horizontal coordinate in tune plots, default = xmin
% xmaxplot: maximum horizontal coordinate in tune plots, default = xmax 
% yminplot: minimum vertical coordinate in tune plots, default = ymin 
% ymaxplot: maximum vertical coordinate in tune plots, default = ymax
% dpminplot: minimum momentum deviation in tune plots, deafult = dpmin
% dpmaxplot: maximum momentum deviation in tune plots, default = dpmax
% x0  : (1xN) array of horizontal values for line plots vs y or dp
% y0  : (1xN) array of horizontal values for line plots vs x or dp
% dp0 : (1XN) array of energy deviation values  for line plots vs x or y
%
% dqx     : horizontal half width of square in tune space for FMA plots,
%           default = 0.001
% dqy     : vertical half width of square in tune space for FMA plots,
%           default = 0.001 
%
% caxrange: color axis range, default [-10 0]
% verbose : defines level of verbose output, default=0, i.e. no output
%
% Optional flags
% plot    : plots tune map
% rate    : chooses rate for diffusion plots
% 
%% Outputs
% Structure with fields
% tunemap.inputs echoes the inputs given to the function 
%   tunemap.inputs.ACHRO
%   tunemap.inputs.plotargs.plottype
%   tunemap.inputs.plotargs.plotmode
%   tunemap.inputs.plotargs.resorder
%   tunemap.inputs.plotargs.qxrange
%   tunemap.inputs.plotargs.qyrange
%   tunemap.inputs.plotargs.xminplot
%   tunemap.inputs.plotargs.xmaxplot
%   tunemap.inputs.plotargs.yminplot
%   tunemap.inputs.plotargs.ymaxplot
%   tunemap.inputs.plotargs.dpminplot
%   tunemap.inputs.plotargs.dpmaxplot
%   tunemap.inputs.plotargs.dqx
%   tunemap.inputs.plotargs.dqy
%   tunemap.inputs.plotargs.x0
%   tunemap.inputs.plotargs.y0
%   tunemap.inputs.plotargs.dp0
%
% tunemap.outputs contains subfields
%   tunemap.outputs.desc          : datetime + input description
%   tunemap.outputs.TMoptions     : echo of input with possible changes
%   tunemap.outputs.nped          : lattice periodicity
%   tunemap.outputs.Qxx           : (1Xnpx) array of horizontal tune values
%   tunemap.outputs.dQxx          : (1Xnpx) array of horizonal tune change values 
%   tunemap.outputs.Qxxfrac       : (1Xnpx) array of fractional horizonal tune values
%   tunemap.outputs.Qyx           : (1Xnpx) array of vertical tune values
%   tunemap.outputs.dQyx          : (1Xnpx) array of vertical tune change values 
%   tunemap.outputs.Qyxfrac       : (1Xnpx) array of fractional vertical tune values
%   tunemap.outputs.Qxy           : (1Xnpy) array of horizontal tune values
%   tunemap.outputs.dQxy          : (1Xnpy) array of horizonal tune change values 
%   tunemap.outputs.Qxyfrac       : (1Xnpy) array of fractional horizonal tune values
%   tunemap.outputs.Qyy           : (1Xnpy) array of vertical tune values
%   tunemap.outputs.dQyy          : (1Xnpy) array of vertical tune change values 
%   tunemap.outputs.Qyyfrac       : (1Xnpy) array of fractional vertical tune values
%   tunemap.outputs.amplx         : (1Xnpx) array of horizontal amplitudes
%   tunemap.outputs.amply         : (1Xnpy) array of vertical amplitudes
%
%   tunemap.outputs.Qxgridxy      : (npx*npyX1) array of horizontal tune values
%   tunemap.outputs.dQxgridxy     : (npx*npyX1) array of horizonal tune change values 
%   tunemap.outputs.Qxgridxyfrac  : (npx*npyX1) array of fractional horizonal tune values
%   tunemap.outputs.Qygridxy      : (npx*npyX1) array of vertical tune values
%   tunemap.outputs.dQygridxy     : (npx*npyX1) array of vertical tune change values 
%   tunemap.outputs.Qygridxyfrac  : (npx*npyx1) array of fractional vertical tune values
%   tunemap.outputs.axgridxy      : (1Xnpx*npy) array of horizontal coordinates
%   tunemap.outputs.aygridxy      : (1Xnpx*npy) array of vertical coordinates
%
%   tunemap.outputs.Qxgridxdp     : (npx*npdX1) array of horizontal tune values
%   tunemap.outputs.dQxgridxdp    : (npx*npdX1) array of horizonal tune change values 
%   tunemap.outputs.Qxgridxdpfrac : (npx*npdX1) array of fractional horizonal tune values
%   tunemap.outputs.Qygridxdp     : (npx*npdX1) array of vertical tune values
%   tunemap.outputs.dQygridxdp    : (npx*npdX1) array of vertical tune change values 
%   tunemap.outputs.Qygridxdpfrac : (npx*npdx1) array of fractional vertical tune values
%   tunemap.outputs.axgridxdp     : (1Xnpx*npd) array of horizontal coordinates
%   tunemap.outputs.dpgridxdp     : (1Xnpx*npd) array of vertical coordinates
%
%   tunemap.outputs.Qxgridydp     : (npy*npdX1) array of horizontal tune values
%   tunemap.outputs.dQxgridydp    : (npy*npdX1) array of horizonal tune change values 
%   tunemap.outputs.Qxgridydpfrac : (npy*npdX1) array of fractional horizonal tune values
%   tunemap.outputs.Qygridydp     : (npy*npdX1) array of vertical tune values
%   tunemap.outputs.dQygridydp    : (npy*npdX1) array of vertical tune change values 
%   tunemap.outputs.Qygridydpfrac : (npy*npdx1) array of fractional vertical tune values
%   tunemap.outputs.aygridydp     : (1Xnpy*npd) array of horizontal coordinates
%   tunemap.outputs.dpgridydp     : (1Xnpy*npd) array of vertical coordinates
%
%   tunemap.outputs.Qxdp        : (1Xndp) array of horizontal tune values
%   tunemap.outputs.dQxdp       : (1Xnpd) array of horizonal tune change values 
%   tunemap.outputs.Qxdpfrac    : (1Xnpd) array of fractional horizonal tune values
%   tunemap.outputs.Qydp        : (1Xndp) array of vertical tune values
%   tunemap.outputs.dQydp       : (1Xnpd) array of vertical tune change values 
%   tunemap.outputs.Qydpfrac    : (1Xnpd) array of fractional vertical tune values
%   tunemap.outputs.dps         : (1Xnpd) array of energy deviations
%
%   tunemap.outputs.Qdifxy    : (npx*npyX1) array of tune diffusions
%   tunemap.outputs.axdifxy   : (npx*npyX1) array of horizontal coordinates for tune diffusion map
%   tunemap.outputs.aydifxy   : (npx*npyX1) array of vertical coordinates for tune diffusion map
%   tunemap.outputs.Qdifxyra  : (npx*npyX1) array of tune diffusion rates 
%
%   tunemap.outputs.Qdifxdp   : (npd*npxX1) array of tune diffusions
%   tunemap.outputs.dpdifxdp  : (npd*npxX1) array of energy deviations for tune diffusion map
%   tunemap.outputs.axdifxdp  : (npd*npxX1) array of vertical coordinates for tune diffusion map
%   tunemap.outputs.Qdifxdpra : (npd*npxX1) array of tune diffusion rates 
%
%   tunemap.outputs.Qdifydp   : (npd*npyX1) array of tune diffusions
%   tunemap.outputs.dpdifydp  : (npx*npyX1) array of energy deviations for tune diffusion map
%   tunemap.outputs.aydifydp  : (npx*npyX1) array of vertical coordinates for tune diffusion map
%   tunemap.outputs.Qdifxdpra : (npx*npyX1) array of tune diffusion rates 
%
%
%% Usage examples
% tunemap = calcTuneMap(ACHRO,[]'plot','desc','Testing...');
% calcTuneMap(ACHRO,[],'nturns',1024,'plot');
% tunemap = calcTuneMap(ACHRO,TMoptions,'plot','xmax',0.007,'mode','x','npx',128);
% tunemap = calcTuneMap(ACHRO,TMoptions,'plot','xmin',-0.007,'xmax',0.004,'mode','x');
% tunemap = calcTuneMap(ACHRO,TMoptions,'xmin',-0.007,'xmax',0.004,'ymin',0.0,'ymax',0.002,'mode','grid');
% tunemap = calcTuneMap(ACHRO,TMoptions,'mode','chro','npd',121,'dpmin',-0.050,'dpmax',0.050);
% tunemap = calcTuneMap(ACHRO,TMoptions,'mode','difxdp','npd',121,'npx', 121, dpmin',-0.050,'dpmax',0.050);
% tunemap = calcTuneMap(ACHRO,TMoptions,'mode','difxy','npx',64+1,'npy',2*64+1,'xmin',-0.005,'xmax',0.005);

%% History
% PFT 2024/04/27: first version, based on calcADTS
% PFT 2024/05/03: included description input and output, added chromatic
%                 tune map
% PFT 2024/05/09: added difusion maps in (x,dp) and (y,dp) planes
% PFT 2024/05/10: added tune calculation on grids in (x,dp) and (y,dp)
%                 planes
% PFT 2024/05/12: added lattice periodicty check
% PFT 2024/05/19: added check that number of turns is larger than 66 for
%                 NAFF
% PFT 2024/07/05: changed input to use a single structure TMoptions
%
%% Input argument parsing
[ACHRO,TMoptions] = getargs(varargin,[],[]);

if (isempty(TMoptions))
   TMoptions.mode='x';
   TMoptions.npx = 11;
   TMoptions.npy = 11;
   TMoptions.npd = 11;
   TMoptions.xmin = -0.006;
   TMoptions.xmax = +0.006;
   TMoptions.ymin = 0.0;
   TMoptions.ymax = 0.004;
   TMoptions.dp   = 0.0;
   TMoptions.dpmin = -0.04;
   TMoptions.dpmax = +0.04;
   TMoptions.nturns = 1024;
   TMoptions.minampx = 30E-6;
   TMoptions.minampy = 30E-6;
   TMoptions.method = 4;
   TMoptions.smooth = false;
end

plotf            = any(strcmpi(varargin,'plot'));
ratef            = any(strcmpi(varargin,'rate'));
verbosef         = getoption(varargin,'verbose',0);
desc             = getoption(varargin,'desc','Tune map calculation');
smooth           = getoption(varargin,'smooth',TMoptions.smooth);
nturns           = getoption(varargin,'nturns',TMoptions.nturns);
mode             = getoption(varargin,'mode',TMoptions.mode);
minampx          = getoption(varargin,'minampx',TMoptions.minampx);
minampy          = getoption(varargin,'minampy',TMoptions.minampy);
xmax             = getoption(varargin,'xmax', TMoptions.xmax);
xmin             = getoption(varargin,'xmin',TMoptions.xmin);
ymax             = getoption(varargin,'ymax',TMoptions.ymax);
ymin             = getoption(varargin,'ymin',TMoptions.ymin);
dp               = getoption(varargin,'dp',TMoptions.dp);
dpmin            = getoption(varargin,'dpmin',TMoptions.dpmin);
dpmax            = getoption(varargin,'dpmax',TMoptions.dpmax);
npx              = getoption(varargin,'npx',TMoptions.npx);
npy              = getoption(varargin,'npy',TMoptions.npy);
npd              = getoption(varargin,'npd',TMoptions.npd);
method           = getoption(varargin,'method',TMoptions.method);

TMoptions.smooth = smooth;
TMoptions.nturns = nturns;
TMoptions.mode   = mode;
TMoptions.minampx= minampx;
TMoptions.minampy= minampy;
TMoptions.xmax   = xmax;
TMoptions.xmin   = xmin;
TMoptions.ymax   = ymax;
TMoptions.ymin   = ymin;
TMoptions.dp     = dp;
TMoptions.dpmin  = dpmin;
TMoptions.dpmax  = dpmax;
TMoptions.npx    = npx;
TMoptions.npy    = npy;
TMoptions.npd    = npd;
TMoptions.method = method;

%
plotmode         = getoption(varargin,'plotmode','abs');
plottype         = getoption(varargin,'plottype','x');
resorder         = getoption(varargin,'resorder',3);
qxrange          = getoption(varargin,'qxrange',[0.0 1.0]);
qyrange          = getoption(varargin,'qyrange',[0.0 1.0]);
xminplot         = getoption(varargin,'xminplot',xmin);
xmaxplot         = getoption(varargin,'xmaxplot',xmax);
yminplot         = getoption(varargin,'yminplot',ymin);
ymaxplot         = getoption(varargin,'ymaxplot',ymax);
dpminplot        = getoption(varargin,'dpminplot',dpmin);
dpmaxplot        = getoption(varargin,'dpmaxplot',dpmax);
x0               = getoption(varargin,'x0','all');
y0               = getoption(varargin,'y0','all');
dp0              = getoption(varargin,'dp0','all');
dqx              = getoption(varargin,'dqx',0.001);
dqy              = getoption(varargin,'dqy',0.001);
caxrange         = getoption(varargin,'caxrange','auto');


%% Preamble
ACHRO  = atdisable_6d(ACHRO);
nped   = atGetRingProperties(ACHRO).Periodicity;  

% If "method" is NAFF or "mode" requires NAFF
% number of turns is a multipe of 6 and the periodicity, otherwise it is
% only a multiple ofthe periodicity
if (method==4 || smooth || strcmpi(mode,'difxy') || strcmpi(mode,'difxdp') || strcmpi(mode,'difydp'))
    ndiv=lcm(nped,6);
else
    ndiv=nped;
end
if (method==4 || smooth || strcmpi(mode,'difxy') || strcmpi(mode,'difxdp') || strcmpi(mode,'difydp'))
    nturns = 2^(log2(nturns));
end
nturns = (fix(nturns/ndiv)+1)*ndiv; 
% if "mathod" is NAFF, minimum number of turns is 66
%
if (method==4 || smooth || strcmpi(mode,'difxy') || strcmpi(mode,'difxdp') || strcmpi(mode,'difydp'))
    nturns = max(nturns,66);
end
Qxx=[];
Qyx=[];
Qxy=[];
Qyy=[];
dQxx=[];
dQyx=[];
dQxy=[];
dQyy=[];
Qxxfrac=[];
Qyxfrac=[];
Qxyfrac=[];
Qyyfrac=[];

amplx=[];
amply=[];

Qxgridxy=[];
Qygridxy=[];
dQxgridxy=[];
dQygridxy=[];
Qxgridxyfrac=[];
Qygridxyfrac=[];

Qxgridxdp=[];
Qygridxdp=[];
dQxgridxdp=[];
dQygridxdp=[];
Qxgridxdpfrac=[];
Qygridxdpfrac=[];

Qxgridydp=[];
Qygridydp=[];
dQxgridydp=[];
dQygridydp=[];
Qxgridydpfrac=[];
Qygridydpfrac=[];

axgridxy=[];
aygridxy=[];

axgridxdp=[];
dpgridxdp=[];

aygridydp=[];
dpgridydp=[];

Qxdp = [];
Qydp = [];
dQxdp = [];
dQydp = [];
Qxdpfrac = [];
Qydpfrac = [];

dps=[];

Qdifxy=[];
axdifxy=[];
aydifxy=[];
Qdifxyra = [];

Qdifxdp=[];
axdifxdp=[];
dpdifxdp=[];
Qdifxdpra=[];

Qdifydp=[];
aydifydp=[];
dpdifydp=[];
Qdifydpra=[];

%% Calculates Tune Map

tstart=tic;
if (verbosef>0)
    fprintf('%s Starting Tune Map calculation ... \n',datetime);
end
switch mode
    case {'x';'X'}
        Rin       = zeros(6,npx);
        amplx     = linspace(xmin,xmax,npx);
        amply     = zeros(1,npx);
        Rin(1,:)  = amplx;
        Rin(5,:)  = dp;
        [~,x0pos] = min(abs(amplx));
        tunes     = calcTune(ACHRO,Rin,'nturns',nturns,...
                    'method',method,'minampx', ...
                    minampx,'minampy',minampy,'fixturns',...
                    'verbose', verbosef-1);
        Qxx       = tunes.outputs.Qx;
        Qyx       = tunes.outputs.Qy;
        dQxx = Qxx - Qxx(x0pos);
        dQyx = Qyx - Qyx(x0pos); 
        Qxxfrac = Qxx-fix(Qxx);
        Qyxfrac = Qyx-fix(Qyx);

    case {'y';'Y'}
        Rin       = zeros(6,npy);
        amply     = linspace(ymin,ymax,npy);
        Rin(3,:)  = amply;
        Rin(5,:)  = dp;
        [~,y0pos] = min(abs(amply));
        tunes = calcTune(ACHRO,Rin,'nturns',nturns,...
                    'method',method,'minampx', minampx, ...
                    'minampy',minampy,'fixturns',...
                    'verbose', verbosef-1);
        Qxy  = tunes.outputs.Qx;
        Qyy  = tunes.outputs.Qy;

        dQxy = Qxy - Qxy(y0pos);
        dQyy = Qyy - Qyy(y0pos);
        Qxyfrac = Qxy-fix(Qxy);
        Qyyfrac = Qyy-fix(Qyy);

    case {'xy';'XY';'xY';'Xy';'td'}
        Rin       = zeros(6,npx+npy);
        amplx     = linspace(xmin,xmax,npx);
        Rin(1,1:npx) = amplx;
        Rin(5,:)  = dp;
        [~,x0pos] = min(abs(amplx));
        amply     = linspace(ymin,ymax,npy);
        Rin(3,npx+1:npx+npy) = amply;
        [~,y0pos] = min(abs(amply));
        tunes = calcTune(ACHRO,Rin,'nturns',nturns,...
                    'method',method,'minampx', minampx, ...
                    'minampy',minampy,'fixturns',...
                    'verbose', verbosef-1);
        Qxx  = tunes.outputs.Qx(1:npx);
        Qyx  = tunes.outputs.Qy(1:npx);
        Qxy  = tunes.outputs.Qx(npx+1:npx+npy);
        Qyy  = tunes.outputs.Qy(npx+1:npx+npy);
     
        dQxx = Qxx - Qxx(x0pos);
        dQyx = Qyx - Qyx(x0pos); 
        Qxxfrac = Qxx-fix(Qxx);
        Qyxfrac = Qyx-fix(Qyx);
     
        dQxy = Qxy - Qxy(y0pos);
        dQyy = Qyy - Qyy(y0pos);
        Qxyfrac = Qxy-fix(Qxy);
        Qyyfrac = Qyy-fix(Qyy);

    case {'gridxy';'GRIDXY'}
        amplx=linspace(xmin,xmax,npx);
        amply=linspace(ymin,ymax,npy);
        [~,x0pos] = min(abs(amplx));
        [~,y0pos] = min(abs(amply));
        x0y0pos = npy*(x0pos-1)+y0pos;
        [amplx_m, amply_m]=meshgrid(amplx,amply);
        axgridxy  = reshape(amplx_m,npx*npy,1);
        aygridxy  = reshape(amply_m,npx*npy,1);
        
        Qxgridxy=nan(npx*npy,1);
        Qygridxy=nan(npx*npy,1);

        if (smooth)
            method=4;
            Rin = zeros(6,npx*npy);
            Rin(1,:) = axgridxy';
            Rin(3,:) = aygridxy';
            Rin(5,:) = dp;
            tunes = calcTune(ACHRO,Rin,'nturns',nturns,'method',method, ...
                             'minampx',minampx,'minampy',minampy,...
                             'fixturns','verbose',verbosef-1);
            Qxgridxy=tunes.outputs.Qx;
            Qygridxy=tunes.outputs.Qy;
        else
            parfor i=1:npx*npy
                tunes = calcTune(ACHRO,[axgridxy(i) 0.0 aygridxy(i) 0.0 dp 0.0]',...
                        'nturns',nturns,'method',method, ...
                        'minampx',minampx,'minampy',minampy,'fixturns',...
                        'verbose',verbosef-1);
                        Qxgridxy(i)=tunes.outputs.Qx;
                        Qygridxy(i)=tunes.outputs.Qy;
            end
        end

        Qxgridxyfrac = Qxgridxy-fix(Qxgridxy);
        Qygridxyfrac = Qygridxy-fix(Qygridxy);
        dQxgridxy   = Qxgridxy - Qxgridxy(x0y0pos);
        dQygridxy   = Qygridxy - Qygridxy(x0y0pos);

    case {'gridxdp';'GRIDXDP'}
        amplx=linspace(xmin,xmax,npx);
        dps=linspace(dpmin,dpmax,npd);
        [dps_m, amplx_m]=meshgrid(dps,amplx);
        dpgridxdp  = reshape(dps_m,npd*npx,1);
        axgridxdp  = reshape(amplx_m,npd*npx,1); 
        [~,x0pos] = min(abs(amplx));
        [~,dp0pos] = min(abs(dps));
        x0dp0pos = npx*(dp0pos-1)+x0pos;
        
        Qxgridxdp =nan(npd*npx,1);
        Qygridxdp =nan(npd*npx,1);
        if (smooth)
            method=4;
            Rin = zeros(6,npd*npx);
            Rin(1,:) = axgridxdp';
            Rin(5,:) = dpgridxdp';
            tunes = calcTune(ACHRO,Rin,'nturns',nturns,'method',method, ...
                             'minampx',minampx,'minampy',minampy,...
                             'fixturns','verbose',verbosef-1);
            Qxgridxdp=tunes.outputs.Qx;
            Qygridxdp=tunes.outputs.Qy;
        else
            parfor i=1:npd*npx
                tunes = calcTune(ACHRO,[axgridxdp(i) 0.0 0.0 0.0 dpgridxdp(i) 0.0]',...
                        'nturns',nturns,'method',method,'minampx', minampx,...
                        'minampy',minampy,'fixturns','verbose',verbosef-1);
                Qxgridxdp(i)=tunes.outputs.Qx;
                Qygridxdp(i)=tunes.outputs.Qy;
            end
        end
        Qxgridxdpfrac = Qxgridxdp-fix(Qxgridxdp);
        Qygridxdpfrac = Qygridxdp-fix(Qygridxdp);
        dQxgridxdp   = Qxgridxdp - Qxgridxdp(x0dp0pos);
        dQygridxdp   = Qygridxdp - Qygridxdp(x0dp0pos);

    case {'gridydp';'GRIDYDP'}
        amply=linspace(ymin,ymax,npy);
        dps=linspace(dpmin,dpmax,npd);
        [dps_m, amply_m]=meshgrid(dps,amply);
        dpgridydp  = reshape(dps_m,npd*npy,1);
        aygridydp  = reshape(amply_m,npd*npy,1);     
        [~,y0pos] = min(abs(amply));
        [~,dp0pos] = min(abs(dps));
        y0dp0pos = npy*(dp0pos-1)+y0pos;

        Qxgridydp =nan(npd*npy,1);
        Qygridydp =nan(npd*npy,1);
        
        if (smooth)
            method=4;
            Rin = zeros(6,npd*npy);
            Rin(3,:) = aygridydp';
            Rin(5,:) = dpgridydp';
            tunes = calcTune(ACHRO,Rin,'nturns',nturns,'method',method, ...
                             'minampx',minampx,'minampy',minampy,...
                             'fixturns','verbose',verbosef-1);
            Qxgridydp=tunes.outputs.Qx;
            Qygridydp=tunes.outputs.Qy;
        else
            parfor i=1:npd*npy
                tunes = calcTune(ACHRO,[0.0 0.0 aygridydp(i) 0.0 dpgridydp(i) 0.0]',...
                        'nturns',nturns,'method',method,'minampx', minampx,...
                        'minampy',minampy,'fixturns','verbose',verbosef-1);
                Qxgridydp(i)=tunes.outputs.Qx;
                Qygridydp(i)=tunes.outputs.Qy;
            end
        end

        Qxgridydpfrac = Qxgridydp-fix(Qxgridydp);
        Qygridydpfrac = Qygridydp-fix(Qygridydp);
        dQxgridydp   = Qxgridydp - Qxgridydp(y0dp0pos);
        dQygridydp   = Qxgridydp - Qxgridydp(y0dp0pos);

    case {'difxy';'DIFXY'}
        method=4;
        amplx=linspace(xmin,xmax,npx);
        amply=linspace(ymin,ymax,npy);
        [amplx_m, amply_m]=meshgrid(amplx,amply);
        axdifxy  = reshape(amplx_m,npx*npy,1);
        aydifxy  = reshape(amply_m,npx*npy,1);     
        [~,x0pos] = min(abs(amplx));
        [~,y0pos] = min(abs(amply));
        x0y0pos = npy*(x0pos-1)+y0pos;

        Qxgridxy =nan(npx*npy,1);
        Qygridxy =nan(npx*npy,1);
        Qxgridxy2=nan(npx*npy,1);
        Qygridxy2=nan(npx*npy,1);

        if (smooth)
            Rin = zeros(6,npx*npy);
            Rin(1,:) = axdifxy';
            Rin(3,:) = aydifxy';
            Rin(5,:) = dp;
            tunes = calcTune(ACHRO,Rin,'nturns',nturns,'method',method, ...
                             'minampx',minampx,'minampy',minampy,...
                             'nsets',2,'fixturns','verbose',verbosef-1);
            Qxgridxy=tunes.outputs.Qx(:,1);
            Qygridxy=tunes.outputs.Qy(:,1);
            Qxgridxy2=tunes.outputs.Qx(:,2);
            Qygridxy2=tunes.outputs.Qy(:,2);
        else
            parfor i=1:npx*npy
                tunes = calcTune(ACHRO,[axdifxy(i) 0.0 aydifxy(i) 0.0 dp 0.0]',...
                            'nturns',nturns,'method',method,'minampx', minampx,...
                            'minampy',minampy,'nsets',2,'fixturns',...
                            'verbose',verbosef-1);
                Qxgridxy(i)=tunes.outputs.Qx(1);
                Qygridxy(i)=tunes.outputs.Qy(1);
                Qxgridxy2(i)=tunes.outputs.Qx(2);
                Qygridxy2(i)=tunes.outputs.Qy(2);
            end
        end
        Qdifxy = log10(sqrt(((Qxgridxy2-Qxgridxy).^2) + ((Qygridxy2-Qygridxy).^2)));
        Qdifxyra = log10(sqrt(((Qxgridxy2-Qxgridxy).^2) + ((Qygridxy2-Qygridxy).^2))/(nturns));
        Qxgridxyfrac = Qxgridxy - fix(Qxgridxy);
        Qygridxyfrac = Qygridxy - fix(Qygridxy);
        dQxgridxy    = Qxgridxy - Qxgridxy(x0y0pos);
        dQygridxy    = Qygridxy - Qygridxy(x0y0pos);
       
    case {'difxdp';'DIFXDP'}
        method=4;
        amplx=linspace(xmin,xmax,npx);
        dps=linspace(dpmin,dpmax,npd);
        [dps_m, amplx_m]=meshgrid(dps,amplx);
        dpdifxdp  = reshape(dps_m,npd*npx,1);
        axdifxdp  = reshape(amplx_m,npd*npx,1); 
        [~,x0pos] = min(abs(amplx));
        [~,dp0pos] = min(abs(dps));
        x0dp0pos = npx*(dp0pos-1)+x0pos;
        
        Qxgridxdp =nan(npd*npx,1);
        Qygridxdp =nan(npd*npx,1);
        Qxgridxdp2=nan(npd*npx,1);
        Qygridxdp2=nan(npd*npx,1);

        if (smooth)
            Rin = zeros(6,npd*npx);
            Rin(1,:) = axdifxdp';
            Rin(5,:) = dpdifxdp';
            tunes = calcTune(ACHRO,Rin,'nturns',nturns,'method',method, ...
                             'minampx',minampx,'minampy',minampy,...
                             'nsets',2,'fixturns','verbose',verbosef-1);
            Qxgridxdp=tunes.outputs.Qx(:,1);
            Qygridxdp=tunes.outputs.Qy(:,1);
            Qxgridxdp2=tunes.outputs.Qx(:,2);
            Qygridxdp2=tunes.outputs.Qy(:,2);
        else
            parfor i=1:npd*npx
                tunes = calcTune(ACHRO,[axdifxdp(i) 0.0 0.0 0.0 dpdifxdp(i) 0.0]',...
                        'nturns',nturns,'method',method,'minampx', minampx,...
                        'minampy',minampy,'nsets',2,'fixturns',...
                        'verbose',verbosef-1);
                Qxgridxdp(i)=tunes.outputs.Qx(1);
                Qygridxdp(i)=tunes.outputs.Qy(1);
                Qxgridxdp2(i)=tunes.outputs.Qx(2);
                Qygridxdp2(i)=tunes.outputs.Qy(2);
            end
        end
        Qdifxdp = log10(sqrt(((Qxgridxdp2-Qxgridxdp).^2) + ((Qygridxdp2-Qygridxdp).^2)));
        Qdifxdpra = log10(sqrt(((Qxgridxdp2-Qxgridxdp).^2) + ((Qygridxdp2-Qygridxdp).^2))/(nturns));
        Qxgridxdpfrac = Qxgridxdp-fix(Qxgridxdp);
        Qygridxdpfrac = Qygridxdp-fix(Qygridxdp);
        dQxgridxdp   = Qxgridxdp - Qxgridxdp(x0dp0pos);
        dQygridxdp   = Qygridxdp - Qygridxdp(x0dp0pos);


    case {'difydp';'DIFYDP'}
        method=4;
        amply=linspace(ymin,ymax,npy);
        dps=linspace(dpmin,dpmax,npd);
        [dps_m, amply_m]=meshgrid(dps,amply);
        dpdifydp  = reshape(dps_m,npd*npy,1);
        aydifydp  = reshape(amply_m,npd*npy,1);     
        [~,y0pos] = min(abs(amply));
        [~,dp0pos] = min(abs(dps));
        y0dp0pos = npy*(dp0pos-1)+y0pos;

        Qxgridydp =nan(npd*npy,1);
        Qygridydp =nan(npd*npy,1);
        Qxgridydp2=nan(npd*npy,1);
        Qygridydp2=nan(npd*npy,1);

        if (smooth)
            Rin = zeros(6,npd*npy);
            Rin(3,:) = aydifydp';
            Rin(5,:) = dpdifydp';
            tunes = calcTune(ACHRO,Rin,'nturns',nturns,'method',method, ...
                             'minampx',minampx,'minampy',minampy,...
                             'nsets',2,'fixturns','verbose',verbosef-1);
            Qxgridydp=tunes.outputs.Qx(:,1);
            Qygridydp=tunes.outputs.Qy(:,1);
            Qxgridydp2=tunes.outputs.Qx(:,2);
            Qygridydp2=tunes.outputs.Qy(:,2);
        else
            parfor i=1:npd*npy
                tunes = calcTune(ACHRO,[0.0 0.0 aydifydp(i) 0.0 dpdifydp(i) 0.0]',...
                        'nturns',nturns,'method',method,'minampx', minampx,...
                    'minampy',minampy,'nsets',2,'fixturns',...
                    'verbose',verbosef-1);
                Qxgridydp(i)=tunes.outputs.Qx(1);
                Qygridydp(i)=tunes.outputs.Qy(1);
                Qxgridydp2(i)=tunes.outputs.Qx(2);
                Qygridydp2(i)=tunes.outputs.Qy(2);
            end
        end
        Qdifydp = log10(sqrt(((Qxgridydp2-Qxgridydp).^2) + ((Qygridydp2-Qygridydp).^2)));
        Qdifydpra = log10(sqrt(((Qxgridydp2-Qxgridydp).^2) + ((Qygridydp2-Qygridydp).^2))/(nturns));
        Qxgridydpfrac = Qxgridydp-fix(Qxgridydp);
        Qygridydpfrac = Qygridydp-fix(Qygridydp);
        dQxgridydp   = Qxgridydp - Qxgridydp(y0dp0pos);
        dQygridydp   = Qygridydp - Qygridydp(y0dp0pos);
    
    case {'chro';'CHRO'}
         Rin       = zeros(6,npd);
         dps       = linspace(dpmin,dpmax,npd);
         Rin(5,:)  = dps;
         [~,dp0pos] = min(abs(dps));
         tunes     = calcTune(ACHRO,Rin,'nturns',nturns,...
                    'method',method,'minampx', minampx,'minampy',minampy,...
                    'fixturns','verbose',verbosef-1);
         Qxdp      = tunes.outputs.Qx;
         dQxdp     = Qxdp - Qxdp(dp0pos);
         Qxdpfrac  = Qxdp - fix(Qxdp);
         Qydp      = tunes.outputs.Qy;
         dQydp     = Qydp - Qydp(dp0pos);
         Qydpfrac  = Qydp - fix(Qydp);

     otherwise
        fprintf('%s Error in calcTuneMap. Unknown calculation mode : %s \n', ...
                   datetime, mode);
        tunemap=[];
        return
end
telapsed = toc(tstart);
%
%% Collects output structure info
tunemap.inputs.ACHRO=ACHRO;
tunemap.inputs.plotargs.plotmode=plotmode;
tunemap.inputs.plotargs.plottype=plottype;
tunemap.inputs.plotargs.resorder=resorder;
tunemap.inputs.plotargs.qxrange=qxrange;
tunemap.inputs.plotargs.qyrange=qyrange;
tunemap.inputs.plotargs.xminplot=xminplot;
tunemap.inputs.plotargs.xmaxplot=xmaxplot;
tunemap.inputs.plotargs.yminplot=yminplot;
tunemap.inputs.plotargs.ymaxplot=ymaxplot;
tunemap.inputs.plotargs.dpminplot=dpminplot;
tunemap.inputs.plotargs.dpmaxplot=dpmaxplot;
tunemap.inputs.plotargs.dqx=dqx;
tunemap.inputs.plotargs.dqy=dqy;
tunemap.inputs.plotargs.caxrange=caxrange;
tunemap.inputs.plotargs.x0=x0;
tunemap.inputs.plotargs.y0=y0;
tunemap.inputs.plotargs.dp0=dp0;
tunemap.inputs.plotargs.ratef=ratef;

tunemap.outputs.desc=strcat(sprintf('%s',datetime),' : ', desc);
tunemap.outputs.nped=nped;
if (TMoptions.nturns~=nturns)
    fprintf('%s Warning: n. turns changed from %5d %5d \n',...
            datetime, TMoptions.nturns,nturns);
    TMoptions.nturns=nturns;
end

tunemap.outputs.TMoptions=TMoptions;
tunemap.outputs.Qxx=Qxx;
tunemap.outputs.Qyx=Qyx;
tunemap.outputs.Qxy=Qxy;
tunemap.outputs.Qyy=Qyy;
tunemap.outputs.dQxx=dQxx;
tunemap.outputs.dQyx=dQyx;
tunemap.outputs.dQxy=dQxy;
tunemap.outputs.dQyy=dQyy;
tunemap.outputs.Qxxfrac=Qxxfrac;
tunemap.outputs.Qyxfrac=Qyxfrac;
tunemap.outputs.Qxyfrac=Qxyfrac;
tunemap.outputs.Qyyfrac=Qyyfrac;

tunemap.outputs.amplx=amplx;
tunemap.outputs.amply=amply;

tunemap.outputs.Qxgridxy=Qxgridxy;
tunemap.outputs.Qygridxy=Qygridxy;
tunemap.outputs.dQxgridxy=dQxgridxy;
tunemap.outputs.dQygridxy=dQygridxy;
tunemap.outputs.Qxgridxyfrac=Qxgridxyfrac;
tunemap.outputs.Qygridxyfrac=Qygridxyfrac;
tunemap.outputs.axgridxy=axgridxy;
tunemap.outputs.aygridxy=aygridxy;

tunemap.outputs.Qxgridxdp=Qxgridxdp;
tunemap.outputs.Qygridxdp=Qygridxdp;
tunemap.outputs.dQxgridxdp=dQxgridxdp;
tunemap.outputs.dQygridxdp=dQygridxdp;
tunemap.outputs.Qxgridxdpfrac=Qxgridxdpfrac;
tunemap.outputs.Qygridxdpfrac=Qygridxdpfrac;
tunemap.outputs.axgridxdp=axgridxdp;
tunemap.outputs.dpgridxdp=dpgridxdp;

tunemap.outputs.Qxgridydp=Qxgridydp;
tunemap.outputs.Qygridydp=Qygridydp;
tunemap.outputs.dQxgridydp=dQxgridydp;
tunemap.outputs.dQygridydp=dQygridydp;
tunemap.outputs.Qxgridydpfrac = Qxgridydpfrac;
tunemap.outputs.Qygridydpfrac = Qygridydpfrac;
tunemap.outputs.aygridydp = aygridydp;
tunemap.outputs.dpgridydp = dpgridydp;

tunemap.outputs.Qxdp     = Qxdp;
tunemap.outputs.dQxdp    = dQxdp;
tunemap.outputs.Qxdpfrac = Qxdpfrac;
tunemap.outputs.Qydp     = Qydp;
tunemap.outputs.dQydp    = dQydp;
tunemap.outputs.Qydpfrac = Qydpfrac;

tunemap.outputs.dps      = dps;

tunemap.outputs.Qdifxy    = Qdifxy;
tunemap.outputs.axdifxy   = axdifxy;
tunemap.outputs.aydifxy   = aydifxy;
tunemap.outputs.Qdifxyra  = Qdifxyra;

tunemap.outputs.Qdifxdp   = Qdifxdp;
tunemap.outputs.axdifxdp  = axdifxdp;
tunemap.outputs.dpdifxdp  = dpdifxdp;
tunemap.outputs.Qdifxdpra = Qdifxdpra;

tunemap.outputs.Qdifydp   = Qdifydp;
tunemap.outputs.aydifydp  = aydifydp;
tunemap.outputs.dpdifydp  = dpdifydp;
tunemap.outputs.Qdifydpra = Qdifydpra;

tunemap.outputs.telapsed  = telapsed;


%% Plots Tune Map
if (plotf)
    if (ratef)
        plotTuneMap(tunemap,'rate');
    else
        plotTuneMap(tunemap);
    end
end



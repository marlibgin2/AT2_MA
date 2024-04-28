function tunemap=calcTuneMap(varargin)
% Calculates and plots betatron tunes maps
% as a function of position or energy.
% Calculates an plots frequency maps and tune diffusion maps
%
% This is a high level wrapper to the function calcTune. 
% Tracking is 5d, i.e. 4d + fixed energy deviation
%   
%% Inputs 
% Mandatory arguments
% RING : AT2 lattice array. 
%
% Optional arguments
% mode: 'X', 'Y', 'XY', 'td', 'grid', 'diff', 'chrom': default is 'X'
%  'x'     calculates only tunes vs horizontal position
%  'y'     calculates only tunes vs vertical position
%  'xy'    calculates both x and y
%  'grid'  calculates tunes on a grid of points in xy plane.
%  'diff'  calculates tune diffusion map
%  'chrom' calculates chromatic tune footprint.
%
% minamplx: minimum absolute value of amplitude in horizontal direction,
%            default=30 microm
% minamply: minimum absolute value of amplitude in horizontal direction,
%            default=30 microm
% xmax: max horizontal amplitude [m], default = 0.005 
% xmin: min horizontal amplitude [m], default = 0.0 
%       
% ymax: max vertical amplitude [m], default = 0.004
% ymin: min vertical amplitude [m], default = 0.0
% 
% dpmax: max energy deviation for chromatic tune footprint,default = +0.03
% dpmin: min energy deviation for chromatic tune footprint,default = -0.03
%
% dp: initial energy deviation, default = 0.0
%
% npx: number of points along horizontal direction; default = 11
% npy: number of points along vertical direction: default = 11
% ndp: number of points along energy axis: default = 11
% nturns : number of turns, default = 128
%
% method 1: Highest peak in fft
%        2: Interpolation on fft results
%        3: Windowing + interpolation (default)
%        4: NAFF
%
% plotmode: 'abs' : plots full tune value (inc. integer part) (default)
%           'rel' : plots tune variations with respect to small amplitude
%                   tunes.
% plottype: 
%   'x'   : ADTS along horizontal axis (default)
%   'y'   : ADTS along vertical axis
%   'xy'  : both ADTS
%   'td'  : same as above, but on a tune diagram
%   'grid': grid of points on a tune diagram 
%   'diff': tune diffusion plot on xy plane
%   'fmap': tune diffusion plot on tune diagram
%   'chro': chromatic tune footprint on a tune diagram

% resorder: resonance order for tune diagram, default = 5
% qxrange=[qxmin qymin]: horizontal plot range in tune diagram,  default =[0 1]
% qyrange=[qymin qymax]: vertical plot range in tune diagram, default= [0 1]
%
% Optional flags
% plot : plots tune map
% verbose: produces verbose output
% 
%% Outputs
% Structure with fields
% tunemap.inputs echoes the inputs given to the function 
%   tunemap.inputs.plane
%   tunemap.inputs.xmin 
%   tunemap.inputs.xmax
%   tunemap.inputs.ymax
%   tunemap.inputs.dp
%   tunemap.inputs.dpmin
%   tunemap.inputs.dpmax
%   tunemap.inputs.npx
%   tunemap.inputs.npy
%   tunemap.inputs.ndp
%   tunemap.inputs.nturns 
%   tunemap.inputs.method 
%   tunemap.inputs.minamplx
%   tunemap.inputs.minamply
%   tunemap.inputs.resorder
%   tunemap.inputs.qxrange
%   tunemap.inputs.qyrange
%
% tunemap.outputs.Qxx: (1xnx) array of horizontal tune values
% tunemap.outputs.Qyx: (1xnx) array of vertical tune values
% tunemap.outputs.Qxy: (1xny) array of horizontal tune values
% tunemap.outputs.Qyy: (1xny) array of vertical tune values
% tunemap.outputs.amplx : (1xn) array of horizontal amplitudes
% tunemap.outputs.amply : (1xn) array of vertical amplitudes
% tunemap.outputs.Qxgrid: (nxxny) array of horizontal tune values
% tunemap.outputs.Qygrid: (nxxny) array of vertical tune values
%
%
%% Usage examples
% tunemap = calcTuneMap(RING,'plot');
% calcTuneMap(RING,'nturns',1024,'plot');
% tunemap = calcTuneMap(RING,'plot','xmax',0.007,'mode','x','npx',128);
% tunemap = calcTuneMap(RING,'plot','xmin',-0.007,'xmax',0.004,'mode','x');
% tunemap = calcTuneMap(RING,'xmin',-0.007,'xmax',0.004,'ymin',0.0,'ymax',0.002,'mode','grid');

%% History
% PFT 2024/04/27: first version, based on calcADTS
%
%% Input argument parsing
[RING] = getargs(varargin,[]);
plotf            = any(strcmpi(varargin,'plot'));
verbosef         = any(strcmpi(varargin,'verbose'));
nturns           = getoption(varargin,'nturns',128);
mode             = getoption(varargin,'mode','x');
minampx          = getoption(varargin,'minampx',30E-6);
minampy          = getoption(varargin,'minampy',30E-6);
xmax             = getoption(varargin,'xmax',0.005);
xmin             = getoption(varargin,'xmin',0.0);
ymax             = getoption(varargin,'ymax',0.004);
ymin             = getoption(varargin,'ymin',0.000);
dp               = getoption(varargin,'dp',0.0);
dpmin            = getoption(varargin,'dpmin',-0.03);
dpmax            = getoption(varargin,'dpmax',+0.03);
npx              = getoption(varargin,'npx',11);
npy              = getoption(varargin,'npy',11);
npd              = getoption(varargin,'npd',11);
method           = getoption(varargin,'method',3);
plotmode         = getoption(varargin,'plotmode','abs');
plottype         = getoption(varargin,'plottype','x');
resorder         = getoption(varargin,'resorder',3);
qxrange          = getoption(varargin,'qxrange',[0.0 1.0]);
qyrange          = getoption(varargin,'qyrange',[0.0 1.0]);

%% Preamble
RING=atdisable_6d(RING);
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

Qxgrid=[];
Qygrid=[];
Qxgridfrac=[];
Qygridfrac=[];

Qxgrid2=[];
Qygrid2=[];
Qxgridfrac2=[];
Qygridfrac2=[];


%% Calculates Tune Map

if (verbosef)
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
        tunes     = calcTune(RING,Rin,'nturns',nturns,...
                    'method',method,'minampx', ...
                    minampx,'minampy',minampy);
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
        tunes = calcTune(RING,Rin,'nturns',nturns,...
                    'method',method,'minampx', minampx, 'minampy',minampy);
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
        tunes = calcTune(RING,Rin,'nturns',nturns,...
                    'method',method,'minampx', minampx, 'minampy',minampy);
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

     case {'grid';'GRID'}
        amplx=linspace(xmin,xmax,npx);
        amply=linspace(ymin,ymax,npy);
        [amplx_m, amply_m]=meshgrid(amplx,amply);
        ax  = reshape(amplx_m,npx*npy,1);
        ay  = reshape(amply_m,npx*npy,1);
        Rin = zeros(6,npx*npy);
        Rin(1,:)=ax;
        Rin(3,:)=ay;
        Rin(5,:)=dp;
        tunes = calcTune(RING,Rin,'nturns',nturns,...
                    'method',method,'minampx', minampx,'minampy',minampy);
     
        Qxgrid = tunes.outputs.Qx;
        Qygrid = tunes.outputs.Qy;

        Qxgridfrac = Qxgrid-fix(Qxgrid);
        Qygridfrac = Qygrid-fix(Qygrid);

        case {'diff';'DIFF'}
             fprintf('%s Error in calcTuneMap. Calculation mode not implemented yet: %s \n', ...
                   datetime, mode);


        case {'chro';'CHRO'}

             fprintf('%s Error in calcTuneMap. Calculation mode not implemented yet: %s \n', ...
                   datetime, mode);
        tunemap=[];

     otherwise
        fprintf('%s Error in calcTuneMap. Unknown calculation mode : %s \n', ...
                   datetime, mode);
        tunemap=[];
        return
end
%
%% Collects output structure info
tunemap.inputs.mode=mode;
tunemap.inputs.dp=dp;
tunemap.inputs.xmin=xmin;
tunemap.inputs.xmax=xmax;
tunemap.inputs.ymin=ymin;
tunemap.inputs.ymax=ymax;
tunemap.inputs.dpmin=dpmin;
tunemap.inputs.dpmax=dpmax;
tunemap.inputs.npx=npx;
tunemap.inputs.npy=npy;
tunemap.inputs.npd=npd;
tunemap.inputs.nturns=nturns; 
tunemap.inputs.method=method;
tunemap.inputs.plotmode=plotmode;
tunemap.inputs.plottype=plottype;
tunemap.inputs.minampx=minampx;
tunemap.inputs.minampy=minampy;
tunemap.inputs.resorder=resorder;
tunemap.inputs.qxrange=qxrange;
tunemap.inputs.qyrange=qyrange;

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

tunemap.outputs.Qxgrid=Qxgrid;
tunemap.outputs.Qygrid=Qygrid;
tunemap.outputs.Qxgridfrac = Qxgridfrac;
tunemap.outputs.Qygridfrac = Qygridfrac;

%% Plots Tune Map
if (plotf)
   plotTuneMap(tunemap,'plottype',plottype,'plotmode',plotmode,'resorder', resorder, ...
                 'qxrange', qxrange, 'qyrange', qyrange);
end



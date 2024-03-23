function adts=calcADTS(varargin)
% Calculates and plots ADTS. This is a high level wrapper to
% the AT2.0 function atnuampl. Tracking is 4d or 6d as defined
% by the input lattice
%   
%% Mandatory input arguments
% RING : AT2 lattice array. 
%
%% Optional input parameters
% plane: 'X', 'Y', 'XY', 'td': default is 'X'
%  'x' calculates only tunes vs horizontal position
%  'y' calculates only tunes vs vertical position
%  'xy' or 'td' calculates both but plots (if requested) are different
%  'xy' will bloot both tujes vs x and y whereas 'td' will plot a tune
%       diagram.
% xmax: max horizontal amplitude [m], default = 0.005 
% xmin: min horizontal amplitude [m], default = 0.0 
%       (note: the function atnuampl replaces zero amplitudes with
%        +30 micrometers)
% np: number of points at which to calculate ADTS default = 10
% ymax: max vertical amplitude [m], default =0.004
% nturns : number of turns, default = 128
% method 1: Highest peak in fft
%        2: Interpolation on fft results
%        3: Windowing + interpolation (default)
%        4: NAFF
% plotmode: 'abs' : plots full tune value (including integer part)
%           'rel' : plots tune variations with respect to small amplitude
%                   tunes.
%% Optional flags
% plot : plots ADTS
% verbose: produces verbose output
% 
%% Output parameters
% Structure with fields
% adts.inputs echoes the inputs given to the function 
%   adts.inputs.plane
%   adts.inputs.xmin 
%   adts.inputs.xmax
%   adts.inputs.ymax
%   adts.inputs.np
%   adts.inputs.nturns 
%   adts.inputs.method 
%
% adts.outputs.
%
% adts.outputs.Qxx: (1xnx) array of horizontal tune values
% adts.outputs.Qyx: (1xnx) array of vertical tune values
% adts.outputs.Qxy: (1xny) array of horizontal tune values
% adts.outputs.Qyy: (1xny) array of vertical tune values
% adts.outputs.amplx : (1xn) array of horizontal amplitudes
% adts.outputs.amply : (1xn) array of vertical amplitudes
%
%
%% Usage examples
% adts = calcADTS(RING,'plot');
% calcADTS(RING,'nturns',1024,'plot');
% adts = calcADTS(RING,'plot','xmax',0.0095839,'plane','x','np',128);

%% History
% 2024/03/17: first version PFT.
% 2024/03/23: added possibility of negative amplitudes, documentation
%             improvements. Changed output argument to structure. Allow
%             plot of changes of tunes, rather then the tunes.
%
%% Input argument parsing
[RING] = getargs(varargin,[]);
plotf            = any(strcmpi(varargin,'plot'));
verbosef         = any(strcmpi(varargin,'verbose'));
nturns           = getoption(varargin,'nturns',128);
plane            = getoption(varargin,'plane','x');
xmax             = getoption(varargin,'xmax',0.005);
xmin             = getoption(varargin,'xmin',0.0);
ymax             = getoption(varargin,'ymax',0.004);
np               = getoption(varargin,'np',10);
method           = getoption(varargin,'method',3);
plotmode         = getoption(varargin,'plotmode','abs');

%% Calculates ADTS

if (verbosef)
    fprintf('%s Starting ADTS calculation ... \n',datetime);
end
switch plane
    case {'x';'X'}
        amplx     = linspace(xmin,xmax,np);
        [~,x0pos] = min(abs(amplx));
        [Qxx,Qyx] = atnuampl(RING,amplx,1,'nturns',nturns,'method',method);
        dQxx = Qxx - Qxx(x0pos);
        dQyx = Qyx - Qyx(x0pos); 
        Qxxfrac = Qxx-fix(Qxx);
        Qyxfrac = Qyx-fix(Qyx);
        Qxy=[];
        Qyy=[];
        dQxy=[];
        dQyy=[];
        Qxyfrac=[];
        Qyyfrac=[];
        amply=[];

    case {'y';'Y'}
        amply = linspace(0,ymax,np);
        y0pos = 1;
        [Qxy,Qyy] = atnuampl(RING,amply,3,'nturns',nturns,'method',method);
        dQxy = Qxy - Qxy(y0pos);
        dQyy = Qyy - Qyy(y0pos);
        Qxx=[];
        Qyx=[];
        Qxyfrac = Qxy-fix(Qxy);
        Qyyfrac = Qyy-fix(Qyy);
        dQxx=[];
        dQyx=[];
        Qxxfrac=[];
        Qyxfrac=[];
        amplx=[];

    case {'xy';'XY';'xY';'Xy';'td'}
        amplx     = linspace(xmin,xmax,np);
        [~,x0pos] = min(abs(amplx));
        amply     = linspace(0,ymax,np);
        y0pos     = 1;
        [Qxx,Qyx] = atnuampl(RING,amplx,1,'nturns',nturns,'method',method);
        dQxx = Qxx - Qxx(x0pos);
        dQyx = Qyx - Qyx(x0pos); 
        Qxxfrac = Qxx-fix(Qxx);
        Qyxfrac = Qyx-fix(Qyx);
        [Qxy,Qyy] = atnuampl(RING,amply,3,'nturns',nturns,'method',method);
        dQxy = Qxy - Qxy(y0pos);
        dQyy = Qyy - Qyy(y0pos);
        Qxyfrac = Qxy-fix(Qxy);
        Qyyfrac = Qyy-fix(Qyy);

     otherwise
        fprintf('%s Error in calcADTS. Unknown plane : %s \n', ...
                   datetime, plane);
        adts=[];
        return
end
%
%% Collects output structure info
adts.inputs.plane=plane;
adts.inputs.xmin=xmin;
adts.inputs.xmax=xmax;
adts.inputs.ymax=ymax;
adts.inputs.np=np;
adts.inputs.nturns=nturns; 
adts.inputs.method=method;
adts.inputs.plotmode=plotmode;

adts.outputs.Qxx=Qxx;
adts.outputs.Qyx=Qyx;
adts.outputs.Qxy=Qxy;
adts.outputs.Qyy=Qyy;
adts.outputs.dQxx=dQxx;
adts.outputs.dQyx=dQyx;
adts.outputs.dQxy=dQxy;
adts.outputs.dQyy=dQyy;
adts.outputs.Qxxfrac=Qxxfrac;
adts.outputs.Qyxfrac=Qyxfrac;
adts.outputs.Qxyfrac=Qxyfrac;
adts.outputs.Qyyfrac=Qyyfrac;

adts.outputs.amplx=amplx;
adts.outputs.amply=amply; 

%% Plots ADTS
if (plotf)
   plotADTS(adts);
end

end

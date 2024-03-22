function [Qxx,Qyx, Qxy, Qyy, amplx, amply]=calcADTS(varargin)
% Calculates and plots ADTS. This is a hgh level wrapper to
% the AT2.0 function atnuampl. Tracking is 4d or 6d depending as defined
% by the input lattice
%   
%% Mandatory input arguments
% RING : AT2 lattice array. 
%
%% Optional input parameters
% plane: 'X', 'Y', or 'XY': default is 'X'
% xmax: max horizontal amplitude [m], default = 0.005 
% xmin: min
% npx: number of points along horizontal axis, default = 10
% ymax: max horizontal amplitude [m], default =0.004
% npy: number of points along horizontal axis, default = 10
% nturns : numbr of turns, default = 1024
% method 1: Highest peak in fft
%        2: Interpolation on fft results
%        3: Windowing + interpolation (default)
%        4: NAFF
%
%% Optional flags
% plot : plots ADTS
% verbose: produces verbose output
% 
%% Output parameters
% Qxx: (1xnx) array of horizontal tune values
% Qyx: (1xnx) array of vertical tune values
% Qxy: (1xny) array of horizontal tune values
% Qyy: (1xny) array of vertical tune values
% amplx : (1xn) array of horizontal amplitudes
% amply : (1xn) array of vertical amplitudes
%% Usage examples
% [Qxx,Qyx,~,~,~,~] = calcADTS(RING,'plot');
% calcADTS(RING,'nturns',1024,'plot');
% [Qxx,Qyx,~,~,amplx,~] = calcADTS(RING,'plot','xmax',0.0095839,'x','npx',200,'method',1);

%% History
% PFT 2024/03/17
%
%% Input argument parsing
[RING] = getargs(varargin,[]);
plotf            = any(strcmpi(varargin,'plot'));
verbosef         = any(strcmpi(varargin,'verbose'));
nturns           = getoption(varargin,'nturns',1024);
plane            = getoption(varargin,'nplane','x');
xmax             = getoption(varargin,'xmax',0.005);
ymax             = getoption(varargin,'ymax',0.004);
npx              = getoption(varargin,'npx',10);
npy              = getoption(varargin,'npy',10);
method           = getoption(varargin,'method',3);

%% preamble

amplx = linspace(0,xmax,npx);
amply = linspace(0,ymax,npy);
if (verbosef)
    fprintf('%s Starting ADTS calculation ',datetime);
end
switch plane
    case {'x';'X'}
        [Qxx,Qyx] = atnuampl(RING,amplx,1,'nturns',nturns,'method',method);
        if (plotf)
            figure; xlim([0,xmax*1000]);
            plot(amplx*1000,Qxx,'-ob');xlabel('X[mm]');ylabel('Qx');
            hold on;
            yyaxis right; 
            plot(amplx*1000,Qyx,'-or');ylabel('Qy');
            grid on;
        end
        Qxy=[];
        Qyy=[];

    case {'y';'Y'}
        [Qxy,Qyy] = atnuampl(RING,amply,3,'nturns',nturns,'method',method);
        if (plotf)
            figure; xlim([0,ymax*1000]); 
            plot(amply*1000,Qxy,'-ob');xlabel('Y[mm]');ylabel('Qx');
            hold on;
            yyaxis right; 
            plot(amply*1000,Qyy,'-ob');ylabel('Qy');
            grid on;
        end
        Qxx=[];
        Qyx=[];
end




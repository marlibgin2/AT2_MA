function adts=calcADTS(varargin)
% Calculates and plots ADTS. This is a high level wrapper to
% the AT2.0 function "atnuampl.m". Tracking is 4d or 6d as defined
% by the input lattice
%   
%% Inputs 
% Mandatory arguments
% RING : AT2 lattice array. 
%
% Optional arguments
% plane: 'X', 'Y', 'XY', 'td', 'grid': default is 'X'
%  'x' calculates only tunes vs horizontal position
%  'y' calculates only tunes vs vertical position
%  'xy' or 'td' calculates both but plots (if requested) are different
%  'grid' calculates tunes on a grid of points in xy plane.
%               'xy' will plot both tunes vs x and y 
%                whereas 'td' and 'grid' will plot a tune diagram.
%
% xmax: max horizontal amplitude [m], default = 0.005 
% xmin: min horizontal amplitude [m], default = 0.0 
%       (note: the function atnuampl replaces zero amplitudes with
%        +30 micrometers)
% ymax: max vertical amplitude [m], default = 0.004
% ymin: min vertical amplitude [m], default = 0.0
% npx: number of points along horizontal direction; default = 11
% npy: number of points along vertical direction: default = 11
% nturns : number of turns, default = 128
% method 1: Highest peak in fft
%        2: Interpolation on fft results
%        3: Windowing + interpolation (default)
%        4: NAFF
%
% plotmode: 'abs' : plots full tune value (including integer part)
%           'rel' : plots tune variations with respect to small amplitude
%                   tunes.
% resorder: resonance order for tune diagram, default = 5
% qxrange=[qxmin qymin]: horizontal plot range in tune diagram,  default =[0 1]
% qyrange=[qymin qymax]: vertical plot range in tune diagram, default= [0 1]
%
% Optional flags
% plot : plots ADTS
% verbose: produces verbose output
% 
%% Outputs
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
% adts.outputs.Qxx: (1xnx) array of horizontal tune values
% adts.outputs.Qyx: (1xnx) array of vertical tune values
% adts.outputs.Qxy: (1xny) array of horizontal tune values
% adts.outputs.Qyy: (1xny) array of vertical tune values
% adts.outputs.amplx : (1xn) array of horizontal amplitudes
% adts.outputs.amply : (1xn) array of vertical amplitudes
% adts.outputs.Qxgrid: (nxxny) array of horizontal tune values
% adts.outputs.Qygrid: (nxxny) array of vertical tune values
%
%
%% Usage examples
% adts = calcADTS(RING,'plot');
% calcADTS(RING,'nturns',1024,'plot');
% adts = calcADTS(RING,'plot','xmax',0.007,'plane','x','npx',128);
% adts = calcADTS(RING,'plot','xmin',-0.007,'xmax',0.004,'plane','x');
% adts = calcADTS(RING,'xmin',-0.007,'xmax',0.004,'ymin',0.0,'ymax',0.002,'plane','grid');

%% History
% 2024/03/17: first version PFT.
% 2024/03/23: added possibility of negative amplitudes, documentation
%             improvements. Changed output argument to structure. Allow
%             plot of changes of tunes, rather then the tunes.
% 2024/03/31: added possibiilty of calculation on a 2d grid of points.
%             added furtehr plot options from plotADTS
% 2024/04/14: adap to nbew version of atnuampl - additional inut parameter
%             is minampl.
%
%% Input argument parsing
[RING] = getargs(varargin,[]);
plotf            = any(strcmpi(varargin,'plot'));
verbosef         = any(strcmpi(varargin,'verbose'));
nturns           = getoption(varargin,'nturns',128);
plane            = getoption(varargin,'plane','x');
minamp           = getoption(varargin,'minamp',30E-6);
xmax             = getoption(varargin,'xmax',0.005);
xmin             = getoption(varargin,'xmin',0.0);
ymax             = getoption(varargin,'ymax',0.004);
ymin             = getoption(varargin,'ymin',0.000);
npx              = getoption(varargin,'npx',11);
npy              = getoption(varargin,'npy',11);
method           = getoption(varargin,'method',3);
plotmode         = getoption(varargin,'plotmode','abs');
resorder         = getoption(varargin,'resorder',5);
qxrange          = getoption(varargin,'qxrange',[0.0 1.0]);
qyrange          = getoption(varargin,'qyrange',[0.0 1.0]);

%% Preamble
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


%% Calculates ADTS

if (verbosef)
    fprintf('%s Starting ADTS calculation ... \n',datetime);
end
switch plane
    case {'x';'X'}
        amplx     = linspace(xmin,xmax,npx);
        [~,x0pos] = min(abs(amplx));
        [Qxx,Qyx] = atnuampl(RING,amplx,1,'nturns',nturns,...
                    'method',method,'minamp', minamp);
        dQxx = Qxx - Qxx(x0pos);
        dQyx = Qyx - Qyx(x0pos); 
        Qxxfrac = Qxx-fix(Qxx);
        Qyxfrac = Qyx-fix(Qyx);

    case {'y';'Y'}
        amply = linspace(ymin,ymax,npy);
        [~,y0pos] = min(abs(amply));
        [Qxy,Qyy] = atnuampl(RING,amply,3,'nturns',nturns,...
                    'method',method,'minamp', minamp);
        dQxy = Qxy - Qxy(y0pos);
        dQyy = Qyy - Qyy(y0pos);
        Qxyfrac = Qxy-fix(Qxy);
        Qyyfrac = Qyy-fix(Qyy);

    case {'xy';'XY';'xY';'Xy';'td'}
        amplx     = linspace(xmin,xmax,npx);
        [~,x0pos] = min(abs(amplx));
        amply     = linspace(ymin,ymax,npy);
        [~,y0pos] = min(abs(amply));
        [Qxx,Qyx] = atnuampl(RING,amplx,1,'nturns',nturns,...
                    'method',method,'minamp', minamp);
        dQxx = Qxx - Qxx(x0pos);
        dQyx = Qyx - Qyx(x0pos); 
        Qxxfrac = Qxx-fix(Qxx);
        Qyxfrac = Qyx-fix(Qyx);
        [Qxy,Qyy] = atnuampl(RING,amply,3,'nturns',nturns,...
                    'method',method,'minamp', minamp);
        dQxy = Qxy - Qxy(y0pos);
        dQyy = Qyy - Qyy(y0pos);
        Qxyfrac = Qxy-fix(Qxy);
        Qyyfrac = Qyy-fix(Qyy);

     case {'grid';'GRID'}
        amplx=linspace(xmin,xmax,npx);amplx(amplx==0)=0.00003;
        amply=linspace(ymin,ymax,npy);amply(amply==0)=0.00003;
        [amplx_m, amply_m]=meshgrid(amplx,amply);
        ax  = reshape(amplx_m,npx*npy,1);
        ay  = reshape(amply_m,npx*npy,1);
        Rin = zeros(6,npx*npy);
        Rin(1,:)=ax;
        Rin(3,:)=ay;
        %
        % from atnuampl
        %
        dp=0.0;
        [~, orbit]=findorbit4(RING, dp);
        [~,nbper]=atenergy(RING);
        [lindata,fractune0]=atlinopt(RING,dp,1:length(RING)+1, 'orbit', orbit);
        tune0=nbper*lindata(end).mu/2/pi;
        offs=[nbper -nbper];
        
        p1=ringpass(RING,Rin,nturns)-orbit(:,ones(1,npx*npy*nturns));
        tunetrack=[findtune(reshape(p1(1,:),npx*npy,nturns)',method);...
        findtune(reshape(p1(3,:),npx*npy,nturns)',method)]';
        [~,k]=min([fractune0-tunetrack(1,:); 1-fractune0-tunetrack(1,:)]);
        np=offs(k);
        offset=round(tune0-np.*tunetrack(1,:));
        tunetrack=np(ones(npx*npy,1),:).*tunetrack + offset(ones(npx*npy,1),:);
     
        Qxgrid = tunetrack(:,1);
        Qygrid = tunetrack(:,2);

        Qxgridfrac = Qxgrid-fix(Qxgrid);
        Qygridfrac = Qygrid-fix(Qygrid);

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
adts.inputs.ymin=ymin;
adts.inputs.ymax=ymax;
adts.inputs.npx=npx;
adts.inputs.npy=npy;
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

adts.outputs.Qxgrid=Qxgrid;
adts.outputs.Qygrid=Qygrid;
adts.outputs.Qxgridfrac = Qxgridfrac;
adts.outputs.Qygridfrac = Qygridfrac;

%% Plots ADTS
if (plotf)
   plotADTS(adts,'plane',plane,'plotmode',plotmode,'resorder', resorder, ...
                 'qxrange', qxrange, 'qyrange', qyrange);
end



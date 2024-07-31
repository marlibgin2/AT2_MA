function DAS=calcDA(varargin)
%
% Calculates and plots Dynamic Aperture. Tracking can be 6d or 4d
% as defined by the input lattice. This is a higherlevel wrapper function
% that in turn calls the lower level function "calcDA_raw"
% 
%% Inputs
% Mandatory arguments
% RING : AT2 lattice array
% DAoptions :Structure containing the following fields:
%               DAmode   = 'grid', 'border', 'smart_in' or 'smart_out'
%               nturns   : number of turns
%               betax0   : horizontal beta for normalization - if NaN, no normalization is done
%               betay0   : vertical beta for normalization - if NaN no normalization is done
%               xmaxdas  : limits of the range in which the DA border is searched
%               xmindas  : limits of the range in which the DA border is searched
%               ymaxdas  : limits of the range in which the DA border is searched
%               dp       : initial dp/p (6d tracking) or fixed dp/p (4d tracking)
%               npd      : number of points along the energy axis for xdp and ydp 
%                          calculation modes, default = 11
%               dpmin    : minimum energy deviation, default = -0.04
%               dpmax    : maximum energy deviation, default = +0.04
%
%
% Parameters for "border", "smart_in" or "smart_out" DA calculation modes
%               r0      : initial guess [m] (only for "border " mode)
%               nang    : number of angular steps
%               z0      : initial longitudinal coordinate (6d tracking). Nan uses synchrnous phase
%               res     : resolution [m]
%               alpha   : da enlargement factor for border search (only for "border" mode)
% 
%
% Parameters for "grid" DA calculation mode
%               XmaxDA  : Horizontal range is -Xmax to Xmax [m]
%               YmaxDA  : Vertical range is 0 to Ymax [m]
%               npdax   : number of grid points in x direction is 2*npdax+1
%               npday   : number of grid points in y direction is  npday+
%
%
%            If DAoptions = [], hard-coded defaults are used.
%            Values of DAoptions fields are overridden if given explicitly 
%            as input in the form ('parameter', value)
%
% Optional arguments
% all fields in DAoptions listed above. 
% desc : descriptive string, default='DA calculation'
% mode : calculation mode, default = 'xy'. Possible values are
%       'xy'  : calculates DA in the xy plane for a given initial energy
%               deviation. Tracking may be 4d + energy or 6d as defined by
%               the input lattice.
%       'xydp': calculate DA in the xp and yp planes. 
%               Tracking is 4d plus fixed
%               energy deviation.
% verbose : defines level of verbose output, default=0, i.e. no output
%
% Optional flags
% plot : plots DA
%
%% Outputs
% DAS: structure with the fields
% DAS.inputs echoes the input parameters with sub-fields
%   DAS.inputs.RING
%   DAS.inputs.DAoptions
%   DAS.inputs.mode
%
% DAS.outputs has sub-fields
%   
%   DAS.outputs.desc       : datetime + input description
%   DAS.outputs.DA         : Dynamic aperture [mm**2]
%   DAS.outputs.DAV        : (DAoptions.nangX2) array of dynamic aperture border coordinates (if DAoptions'border' DAmode) 
%                            or ( (2*DAoptions.npdax+1)X (DAoptions.npday+1) array of booleans indicating particle loss (if 'grid' DAmode)
%   DAS.outputs.DAoptions  : Structure with options used in the calculation (may differ
%                            from the input DAoptions if specific fields were 
%                            overwritten through optional parameters
%   DAS.outputs.DAXp       : (1 X npd) array of positive X DA [m]
%   DAS.outputs.DAYp       : (1 X npd) array of positive Y DA [m]
%   DAS.outputs.DAXm       : (1 X npd) array of negative X DA [m]
%   DAS.outputs.dps        : (1 X npd) array of energy deviations (mode 'xydp')
%
%   DAS.outputs.telapsed    : elapsed calculation time (sec)
%
%% Usage examples
% DAS = calcDA(RING,DAoptions,'plot');
% DAS = calcDA(RING,[],'desc','Standard Lattice','nturns',1024,'DAmode','grid');
% DAS = calcDA(RING,[],'nturns',1024,'DAmode','grid');
% calcDA(RING,[],'nturns',1024,'DAmode','grid');

%% History
% PFT 2024/03/09
% PFT 2024/05/01 : updated documentation,changed output parameter 
%                  structure and separated calculation from plotting
% PFT 2024/05/13 : added xydp calculation mode
% PFT 2024/05/25 : changed input of xdp and ydp input parameters, which
%                  are now part of the DAoptions structure
%
% PFT 2024/07/28 : adapted to incude DAmode='smart_in'
% PFT 2024/07/30 : added handling of nan as input value for DAoptions.nturns
%
%% Input argument parsing
[RING,DAoptions] = getargs(varargin,[],[]);
DAS.inputs.RING=RING;
DAS.inputs.DAoptions=DAoptions;

if (isempty(DAoptions))
    DAoptions.dp=0.0;
    DAoptions.z0=nan;
    DAoptions.DAmode='grid';
    DAoptions.nturns=nan;
    DAoptions.betax0=nan;
    DAoptions.betay0=nan;
    DAoptions.xmaxdas=0.015;
    DAoptions.xmindas =-0.0150;
    DAoptions.ymaxdas = 0.007;
    DAoptions.npd = 11;
    DAoptions.dpmin = -0.04;
    DAoptions.dpmax = +0.04;
    DAoptions.XmaxDA = 0.015;
    DAoptions.YmaxDA = 0.007;
    DAoptions.npdax = 64;
    DAoptions.npday = 64;
    DAoptions.r0 = 0.020;
    DAoptions.nang = 20;
    DAoptions.res = 5E-4; 
    DAoptions.alpha = 1.1;
    DAoptions.dx=nan;
    DAoptions.dy=nan;
    DAoptions.dxdy=nan;
    DAoptions.npDA = (2*DAoptions.npdax+1)*(DAoptions.npday+1); %total number of grid points
    DAoptions.X0da = zeros(DAoptions.npDA,1);  % horizontal coordinates of grid points [m]
    DAoptions.Y0da = zeros(DAoptions.npDA,1);  % vertical coordinates of grid points [m]
end

mode             = getoption(varargin,'mode','xy');
desc             = getoption(varargin,'desc','DA calculation');
plotf            = any(strcmpi(varargin,'plot'));
verboselevel     = getoption(varargin,'verbose',0);
dp               = getoption(varargin,'dp',DAoptions.dp);
z0               = getoption(varargin,'z0',DAoptions.z0);
DAmode           = getoption(varargin,'DAmode',DAoptions.DAmode);
nturns           = getoption(varargin,'nturns',DAoptions.nturns);
betax0           = getoption(varargin,'betax0',DAoptions.betax0);
betay0           = getoption(varargin,'betay0',DAoptions.betay0);
xmaxdas          = getoption(varargin,'xmaxdas',DAoptions.xmaxdas);
xmindas          = getoption(varargin,'xmindas',DAoptions.xmindas);
ymaxdas          = getoption(varargin,'ymaxdas',DAoptions.ymaxdas);
XmaxDA           = getoption(varargin,'XmaxDA',DAoptions.XmaxDA);
YmaxDA           = getoption(varargin,'YmaxDA',DAoptions.YmaxDA);
npdax            = getoption(varargin,'npdax',DAoptions.npdax);
npday            = getoption(varargin,'npday',DAoptions.npday);
r0               = getoption(varargin,'r0', DAoptions.r0);
nang             = getoption(varargin,'nang',DAoptions.nang);
res              = getoption(varargin,'res',DAoptions.res); 
alpha            = getoption(varargin,'alpha',DAoptions.alpha);
npd              = getoption(varargin,'npd',DAoptions.npd);
dpmin            = getoption(varargin,'dpmin',DAoptions.dpmin);
dpmax            = getoption(varargin,'dpmax',DAoptions.dpmax);

dps=linspace(dpmin,dpmax,npd);

DAS.inputs.mode=mode;

DAoptions.dp=dp;
DAoptions.z0=z0;
DAoptions.DAmode=DAmode;
DAoptions.nturns=nturns;
DAoptions.betax0=betax0;
DAoptions.betay0=betay0;
DAoptions.xmaxas=xmaxdas;
DAoptions.xmindas=xmindas;
DAoptions.ymaxas=ymaxdas;
DAoptions.npd=npd;
DAoptions.dpmin=dpmin;
DAoptions.dpmax=dpmax;
DAoptions.XmaxDA=XmaxDA;
DAoptions.YmaxDA=YmaxDA;
DAoptions.npdax=npdax;
DAoptions.npday=npday;
DAoptions.r0=r0;
DAoptions.nang=nang;
DAoptions.res=res;
DAoptions.alpha=alpha;

% Parameters for dynamic aperture calculation
%

%% Checks backward compatibilty
if (~isfield(DAoptions,'xmaxdas'))
    DAoptions.xmaxdas = 0.007;
end

if (~isfield(DAoptions,'xmindas'))
    DAoptions.xmindas = -0.015;
end

if (~isfield(DAoptions,'ymaxdas'))
    DAoptions.ymaxdas = 0.003;
end
   
if(~(isfield(DAoptions,'XmaxDA')))
    DAoptions.XmaxDA = 0.015;
end

if(~(isfield(DAoptions,'YmaxDA')))
    DAoptions.YmaxDA     = 0.004;
end

XmaxDA = DAoptions.XmaxDA;
YmaxDA = DAoptions.YmaxDA;

if(~(isfield(DAoptions,'DAmode')))
    DAoptions.DAmode = 'border';
end
DAmode = DAoptions.DAmode;

if (strcmp(DAmode,'grid'))
%
%% Recalculates X0da and Y0da in case the data in DAoptions 
% is not consistent (e.g. is not the same as in a previous MOGA run)
%
    npdax      = DAoptions.npdax;   % number of grid points in x direction is 2*npdax+1
    npday      = DAoptions.npday;   % number of grid points in y direction is  npday+1
    npDA       = DAoptions.npDA;    % total numbr of grid points
    dx = XmaxDA/npdax; % grid stepsize in x [m]
    dy = YmaxDA/npday; % grid stepsize in y [m]
    dxdy = dx*dy;
    X0da = zeros(npDA,1);  % horizontal coordinates of grid points [m]
    Y0da = zeros(npDA,1);  % vertical coordinates of grid points [m]
    k= 1;
    for i=0:npday 
        for j=0:2*npdax
        X0da(k) = -XmaxDA+dx*j;
        Y0da(k) =  dy*i;
        k=k+1;
        end
    end
    DAoptions.X0da=X0da;
    DAoptions.Y0da=Y0da;
    DAoptions.dxdy=dxdy;
end
%% preamble
PC=load('PC.mat');      %to prevent matlab from complaining about variable name being the same as script name.
PhysConst = PC.PC;      %Load physical constants

DAXp=zeros(1,npd);
DAXm=zeros(1,npd);
DAYp=zeros(1,npd);

%% Calculates DA
tstart=tic;
if (verboselevel>0)
    fprintf('*** \n');
    fprintf('%s Starting DA calculation \n', datetime);
end
try
   rpara = atsummary(RING);
   etax = rpara.etax;
   % for 6d tracking and if DAoptions.z0 is not given
   % make sure the initial longitudinal coordinate corresponds 
   % to the synchronous particle phase
   if (isnan(z0))
       if (check_6d(RING))
        z0 = PhysConst.c*(rpara.syncphase-pi)/(2*pi*rpara.revFreq*rpara.harmon); %if z0 not given, choose the synchronous phase
       else
        z0=0.0;
       end
       DAoptions.z0=z0;
   end
   %
   % if input number of turns is nan, and lattice is 6d set it to
   %  1.2*synchrotron period
   if (isnan(nturns))
       if (check_6d(RING))
            DAoptions.nturns = round(1.2/rpara.synctune);
       else
            DAoptions.nturns = 1024;
       end
       DAoptions.nturns=nturns;
   end

   switch mode
       case {'xy';'XY'}
            
            [DA,DAV] = calcDA_raw(RING,DAoptions,etax,rpara.beta0(1),rpara.beta0(2));

       case {'xydp';'XYDP'}
            DA=nan;
            if (strcmpi(DAoptions.DAmode,'grid'))
                DAoptions.DAmode = 'smart_in';
            end
            DAoptions.nang = 2;
            for i=1:npd
                if (verboselevel>0)
                    fprintf('%s dp = %4.2f %% \n', datetime, dps(i)*100);
                end
                DAoptions.dp=dps(i);
                [~,DAV]=calcDA_raw(RING,DAoptions,etax,rpara.beta0(1),rpara.beta0(2));
                DAXp(i)=DAV(1,1);
                DAYp(i)=DAV(2,2);
                DAXm(i)=DAV(3,1);
            end

       otherwise
           fprintf('%s Error in calcDA, unknown mode %s \n',datetime, mode);
           DA=NaN;
           DAV=[];
   end

   telapsed=toc(tstart);

catch ME
     fprintf('%s Error in calcDA \n', datetime);
     fprintf('Error message was:%s \n',ME.message);
end

%% Collects data for output structure

   DAS.outputs.desc=strcat(sprintf('%s',datetime),' : ', desc);
   DAS.outputs.DA=DA;
   DAS.outputs.DAV=DAV;
   DAS.outputs.DAoptions=DAoptions;
   DAS.outputs.DAXp=DAXp;
   DAS.outputs.DAYp=DAYp;
   DAS.outputs.DAXm=DAXm;
   DAS.outputs.dps=dps;
   DAS.outputs.telapsed=telapsed;

if (plotf&&not(isempty(DAS)))
    plotDA(DAS);
end


if(verboselevel>0)
    fprintf('%s calcDA: DA calculation complete \n', datetime);
end

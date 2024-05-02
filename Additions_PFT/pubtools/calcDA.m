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
%               DAmode   = 'grid' or 'border'
%               nturns   : number of turns
%               betax0   : horizontal beta for normalization - if NaN, no normalization is done
%               betay0   : vertical beta for normalization - if NaN no normalization is done
%               xmaxdas  : limits of the range in which the DA border is searched
%               xmindas  : limits of the range in which the DA border is searched
%               ymaxdas  : limits of the range in which the DA border is searched
%
% Parameters for "border" DA calculation mode
%               r0      : initial guess [m]
%               nang    : number of angular steps
%               dp      : initial dp/p (6d tracking) or fixed dp/p (4d tracking)
%               z0      : initial longitudinal coordinate (6d tracking). Nan uses synchrnous phase
%               res     : resolution [m]
%               alpha   : da enlargement factor for border search
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
% all fields in DAoptions listed above. If given
%
% Optional flags
% plot : plots DA;
% verbose: produces verbose output
%
%% Outputs
% DAS: structure with the fields
% DAS.inputs echoes the input parameters with fields
%   DAS.inputs.RING
%   DAS.inputs.DAoptions
%
% DAS.outputs has fiels
%   DA       : Dynamic aperture [mm**2]
%   DAV      : Vector of dynamic aperture border coordinates (if 'border' mode) 
%              or vector of booleans indicating particle loss (if 'grid' mode)
%   DAoptions: Structure with options used in the calculation (may differ
%              from the input DAoptions if specific fields were 
%              overwritten through optional parameters
%
%% Usage examples
% DAS = calcDA(RING,DAoptions,'plot');
% DAS = calcDA(RING,[],'nturns',1024,'DAmode','grid');
% calcDA(RING,[],'nturns',1024,'DAmode','grid');

%% History
% PFT 2024/03/09
% PFT 2024/05/01 : updated documentation,changed output parameter 
%                  structure and separated calculation from plotting
%
%% Input argument parsing
[RING,DAoptions] = getargs(varargin,[],[]);
DAS.inputs.DAoptions=DAoptions;
DAS.inputs.RING=RING;
if (isempty(DAoptions))
    DAoptions.dp=0.0;
    DAoptions.z0=nan;
    DAoptions.DAmode='grid';
    DAoptions.nturns=1024;
    DAoptions.betax0=nan;
    DAoptions.betay0=nan;
    DAoptions.xmaxdas=0.007;
    DAoptions.xmindas =-0.0150;
    DAoptions.ymaxdas = 0.006;
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

plotf            = any(strcmpi(varargin,'plot'));
verbosef         = any(strcmpi(varargin,'verbose'));
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

DAoptions.dp=dp;
DAoptions.z0=z0;
DAoptions.DAmode=DAmode;
DAoptions.nturns=nturns;
DAoptions.betax0=betax0;
DAoptions.betay0=betay0;
DAoptions.xmaxas=xmaxdas;
DAoptions.xmindas=xmindas;
DAoptions.ymaxas=ymaxdas;
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
PC=load('PC.mat');      %to prevent matlab from complaining about variable name being the same as script name.
PhysConst = PC.PC;      %Load physical constants

%% Calculates DA
if (verbosef)
    tic;
    fprintf('*** \n');
    fprintf('%s Starting DA calculation \n', datetime);
end
try
   rpara = atsummary(RING);
   etax = rpara.etax;
   %
   % for 6d tracking and if DAoptions.z0 is not given
   % make sure the initial logiudinal coordinate corresponds 
   % to the synchronous particle phase
   if (isnan(z0))
       if (check_6d(RING))
        z0 = PhysConst.c*(rpara.syncphase-pi)/(2*pi*rpara.revFreq*rpara.harmon); %if z0 not given choose the synchronous phase
       else
        z0=0.0;
       end
       DAoptions.z0=z0;
   end
   [DA,DAV] = calcDA_raw(RING,DAoptions,etax,rpara.beta0(1),rpara.beta0(2));
   %% Collects data for output structure
    DAS.outputs.DA=DA;
    DAS.outputs.DAV=DAV;
    DAS.outputs.DAoptions=DAoptions;
catch ME
     fprintf('%s Error in calcDA \n', datetime);
     fprintf('Error message was:%s \n',ME.message);
     DAS=[];
end

if (plotf&&not(isempty(DAS)))
    plotDA(DAS);
end

if(verbosef)
    fprintf('DA calculation complete \n');
    toc;
end

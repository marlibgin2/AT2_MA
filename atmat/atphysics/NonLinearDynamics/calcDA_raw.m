function [DA,DAV] = calcDA_raw(RING,DAoptions,etax,betax,betay)
% Calculates Dynamic Aperture by tracking. Tracking can be 6d or 4d
% as defined by the input lattice. This is a lower level function
% that can be called by the higher level wrapper function "calcDA"  
% or by optimization functions in MOGA/SOGA.
%
%% Inputs
% Mandatory arguments
% RING: AT2 lattice array
% DAoptions: Structure containing the following fields:
%               DAmode   = 'grid', 'border' or 'smart'
%               nturns   : number of turns
%               betax0   : horizontal beta for normalization - if NaN, no normalization is done
%               betay0   : vertical beta for normalization - if NaN no normalization is done
%               xmaxdas  : limits of the range in which the DA border is searched
%               xmindas  : limits of the range in which the DA border is searched
%               ymaxdas  : limits of the range in which the DA border is searched
%
%           Parameters for "border" and "smart" DA calculation modes
%               r0      : initial guess [m]
%               nang    : number of angular steps
%               dp      : initial dp/p (6d tracking) or fixed dp/p (4d tracking)
%               z0      : initial longitudinal coordinate (6d tracking). Nan uses synchrnous phase
%               res     : resolution [m]
%               alpha   : da enlargement factor for border search
%
%           Parameters for "grid" DA calculation mode
%               XmaxDA  : Horizontal range is -Xmax to Xmax [m]
%               YmaxDA  : Vertical range is 0 to Ymax [m]
%               npdax   : number of grid points in x direction is 2*npdax+1
%               npday   : number of grid points in y direction is  npday+
%
% etax  : dispersion at the tracking point [m]
% betax : horizontal beta at the tracking point [m]. Used for scaling. If
%         nan, no scaling is done
% betay : vertical beta at the tracking point [m]. USed for scaling. If
%         nan, no scaling is done
%
%% Outputs
% DA: Dynamic aperture [mm**2]
% DAV : vector of dynamic aperture border coordinates (if 'border' mode) 
%       or vector of booleans indicating particle loss (if 'grid' mode)
%
%% Usage examples
% [DA,~]    = calcDA_raw(RING,DAoptions,0.0,9.0,2.0);
% [DA, DAV] = calcDA_raw(RING,DAoptions,0.0,9.0,2.0);
% [DA, DAV] = calcDA_raw(RING,DAoptions,0.0,nan,nan);

%% History 
% PFT 2024/03/09 
% PFT 2024/06/03: added 'smart' DA calculation mode
% PFT 2024/06/13: removed "smart" mode and added "smart_in" and "smart_out"
%                 modes
% PFT 2024/06/28: improved documentation
%
%% Collects data from DAoptions structure
betax0     = DAoptions.betax0;  % hor beta for normalization - if NaN no normalizatinis done
betay0     = DAoptions.betay0;  % ver beta for normalization - if NaN no normalization is done
DAmode     = DAoptions.DAmode;  % DAmode  = 'border', 'grid', 'smart_in' or 'smart_out;
xmaxdas    = DAoptions.xmaxdas;
xmindas    = DAoptions.xmindas;
ymaxdas    = DAoptions.ymaxdas;

% Parameters for "border" DA calculation mode
r0         = DAoptions.r0;      % initial guess [m];
nangs      = DAoptions.nang;    % number of angular steps
dang       = pi/(nangs-1);      % angular step size [rad}
nturns     = DAoptions.nturns;  % number of turns
dp         = DAoptions.dp;      % initial dp/p (6d tracking) or fixed dp/p (4d tracking)
if (check_6d(RING))
    z0     = DAoptions.z0;      % initial longitudinal coordinate (6d tracking)
else
    z0     = 0.0;
end
res        = DAoptions.res;     % resolution [m]
DAalpha    = DAoptions.alpha;   % DA radius enlargement factor

% Parameters for "grid" dynamic aperture calculation
% Warning the vectors X0da, Y0da must be consistent with the other
% parameters. This is guaranteed when using DAoptions fom LatticeOptData 
% after running m4U.m or when using the calcDA functions. ExMOGA
% also checks for that.
%

dxdy       = DAoptions.dxdy;    % grid cell area [m**2]
X0da       = DAoptions.X0da;    % horizontal coordinates of grid points [m]
Y0da       = DAoptions.Y0da;    % vertical coordinates of grid points [m]

%% Calculates Dynamic Aperture

try
  switch DAmode
     case 'border'
       DAV = calcDA_border(RING, r0, nangs, nturns, dp, z0, etax, res, DAalpha, xmaxdas, xmindas, ymaxdas);
       if (not(isnan(betax0))&&not(isnan(betay0))&&not(isnan(betax))&&not(isnan(betay))) 
           DAVN(:,1)=DAV(:,1)*sqrt(betax0/betax);
           DAVN(:,2)=DAV(:,2)*sqrt(betay0/betay);
       else
           DAVN = DAV;
       end
       RDA = DAVN.*DAVN*1E6; % converts to mm**2
       RA  = RDA(:,1)+RDA(:,2);
       DA  = sum(RA)*dang/2;

     case 'grid'
        DAV = calcDA_grid(RING, X0da, Y0da, nturns, dp, z0, etax, xmaxdas, xmindas, ymaxdas);
        DA = sum(DAV)*dxdy;
        if (not(isnan(betax0))&&not(isnan(betay0))&&not(isnan(betax))&&not(isnan(betay))) 
           DA=DA*sqrt(betax0/betax)*sqrt(betay0/betay);
        end 
        DA=DA*1e6; % converts to mm**2

     case 'smart_in'
        DAV = calcDA_smart_in(RING, nangs, nturns, dp, z0, res, xmaxdas, xmindas, ymaxdas);
        if (not(isnan(betax0))&&not(isnan(betay0))&&not(isnan(betax))&&not(isnan(betay))) 
           DAVN(:,1)=DAV(:,1)*sqrt(betax0/betax);
           DAVN(:,2)=DAV(:,2)*sqrt(betay0/betay);
        else
           DAVN = DAV;
        end
        RDA = DAVN.*DAVN*1E6; % converts to mm**2
        RA  = RDA(:,1)+RDA(:,2);
        DA  = sum(RA)*dang/2; 

     case 'smart_out'
        DAV = calcDA_smart_out(RING, nangs, nturns, dp, z0, res, xmaxdas, xmindas, ymaxdas);
        if (not(isnan(betax0))&&not(isnan(betay0))&&not(isnan(betax))&&not(isnan(betay))) 
           DAVN(:,1)=DAV(:,1)*sqrt(betax0/betax);
           DAVN(:,2)=DAV(:,2)*sqrt(betay0/betay);
        else
           DAVN = DAV;
        end
        RDA = DAVN.*DAVN*1E6; % converts to mm**2
        RA  = RDA(:,1)+RDA(:,2);
        DA  = sum(RA)*dang/2; 
        
     otherwise
        fprintf('%s Error in calcDA_raw: uknown DA calculation mode %s \n', datetime, DAmode);
        DA=0.0;
  end

catch ME
    fprintf('%s Error in calcDA_raw:  \n', datetime);
    fprintf('Error message was:%s \n',ME.message);
    error_line = ME.stack(1).line;
    file = ME.stack(1).file;
    fnct = ME.stack(1).name;
    fprintf('at line number %3d \n', error_line);
    fprintf('file %s \n', file);
    fprintf('function %s \n', fnct);
    DA=0.0;
end
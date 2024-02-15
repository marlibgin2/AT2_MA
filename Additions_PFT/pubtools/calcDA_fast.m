function [DynAp,DAV] = calcDA_fast(LAT,DAoptions,etax,betax,betay)
% Calculates Dynamic Aperture by tracking of lattice "LAT"
% input paramtes are defined in the DAoptions structure.
% 
%
% Parameters for dynamic aperture calculation
%
betax0     = DAoptions.betax0;  % hor beta for normalization - if NaN no normalizatinis done
betay0     = DAoptions.betay0;  % ver beta for normalization - if NaN no normalization is done
DAmode     = DAoptions.DAmode;  % DAmode  = 'border' or 'grid';
xmaxdas    = DAoptions.xmaxdas;
xmindas    = DAoptions.xmindas;
ymaxdas    = DAoptions.ymaxdas;

% Parameters for "border" DA calculation mode
r0         = DAoptions.r0;      % initial guess [m];
nangs      = DAoptions.nang;    % number of angular steps
dang       = pi/(nangs-1);                     % angular step size [rad}
nturns     = DAoptions.nturns;  % number of turns
dp         = DAoptions.dp;      % initial dp/p (6d tracking) or fixed dp/p (4d tracking)
z0         = DAoptions.z0;      % initial longitudinal coordinate (6d tracking)
res        = DAoptions.res;     % resolution [m]
DAalpha    = DAoptions.alpha;   % DA radius enlargement factor

% Parameters for "grid" dynamic aperture calculation
% Warning the vectors X0da, Y0da must be consistent with the other
% parameters. This is guaranteed when using DAoptions fom LatticeOptData 
% after running m4U.m or when using the CalcPlotDA functions. ExMOGA
% also checks for that.
%

dxdy       = DAoptions.dxdy;    % grid cell area [m**2]
X0da       = DAoptions.X0da;    % horizontal coordinates of grid points [m]
Y0da       = DAoptions.Y0da;    % vertical coordinates of grid points [m]

try
  switch DAmode
     case 'border'
       DAV = modelDA_sim_par(LAT, r0, nangs, nturns, dp, z0, etax, res, DAalpha, xmaxdas, xmindas, ymaxdas);
       if (not(isnan(betax0))&&not(isnan(betay0))&&not(isnan(betax))&&not(isnan(betay))) 
           DAVN(:,1)=DAV(:,1)*sqrt(betax0/betax);
           DAVN(:,2)=DAV(:,2)*sqrt(betay0/betay);
       else
           DAVN = DAV;
       end
       RDA = DAVN.*DAVN*1E6; % converts to mm**2
       RA  = RDA(:,1)+RDA(:,2);
       DynAp  = sum(RA)*dang/2;

     case 'grid'
        DAV = calcDA_grid(LAT, X0da, Y0da, nturns, dp, z0, etax, xmaxdas, xmindas, ymaxdas);
        DynAp = sum(DAV)*dxdy;
        if (not(isnan(betax0))&&not(isnan(betay0))&&not(isnan(betax))&&not(isnan(betay))) 
           DynAp=DynAp*sqrt(betax0/betax)*sqrt(betay0/betay);
        end 
        DynAp=DynAp*1e6; % converts to mm**2
     otherwise
        fprintf('Uknown DA calculation mode %s \n', DAmode);
        DynAp=0.0;
  end

catch ME
    fprintf('Error calculating Dynamic Aperture \n');
    fprintf('Error message was:%s \n',ME.message);
    DynAp=0.0;
end
function DAV = calcDA_grid(RING, X0, Y0, nturns, dp, z0, etax, xmax, xmin, ymax)
%
% Calculates Dynamic Aperture by tracking particles on a 2D grid of points.
% Tracking can be 6d or 4d as defined by the input lattice. 
% This is a low level function called by calcDA_raw, which in turn is
% called by the higher level wrapper function "calcDA" or by
% optimzation functions used in MOGA/SOGA
%
%% Mandatory input arguments
% X0: (nx1) vector of initial horizontal coordinates [m]
% Y0: (nx1) vector of initial vertical coordinates [m]
%
% nturns: Numbers of turns to track   
% dp: initial momentum or fixed momentum (if tracking is 4d)
% z0: initial longitudinal position [m]
% etax: dispersion function at the tracking point [m]
% xmax: limits of the region where to look for the DA [m]
% xmin: limits of the region where to look for the DA [m]
% ymax: limits of the region where to look for the DA [m]
%
%% Output parameters
% DAV: (nx1) vector containing 
%                   1 if particle is not lost in nturns
%                   0 if particle is lost in n turns
%
%% Usage example
% X0=[0.005 0.0]'; Y0=[0.0 0.0]';
% nturns=1024; dp=0.0; z0=0.0; 
% etax=0.0;
% xmax=.0015;xmin=-0.015;ymax=0.007;
% DAV = calcDA_grid(RING,X0,Y0,nturns,dp,z0,etax,xmax,xmin,ymax)
%

%% History
% 2023/11/06 Written by Pedro F. Tavares
% 2024/02/09 general updates, documetnation and standardization of names
% 2024/03/10: eliminated the chromatic orbit check. the etax input is no
%             longer required and will also be eliminated

%% Calculates Dynamic Aperture
np   = size(X0,1);
DAV  = zeros(np,1);
% Evaluate the Chromatic orbit
% twiss=  gettwiss(THERING, 0.0);

%x0=twiss.etax(1)*dp;
%x = etax*dp;
%Check that the chromatic orbit is stable
%[~, loss]=ringpass(RING,[x 0.0 0 0.0 dp z0]',nturns);
%if (loss)
%    disp('The chromatic closed orbit is not stable. No DA found');
%    DAV =[];
%else
    parfor i=1:np
        x = X0(i);
        y = Y0(i);
        if (x>xmax || x< xmin || y >ymax)
            loss=true;
        else
            [~, loss]=ringpass(RING,[x 0.0 y 0.0 dp z0]',nturns);
        end
        DAV(i)=not(loss); %
    end
%end

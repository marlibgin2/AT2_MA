function DAV = calcDA_border( RING, r0, nsteps, nturns, dp, z0, etax, res, alpha, xmax, xmin, ymax)
%
% Calculates Dynamic Aperture by tracking particles alogm radial lines
% on the XY plane until they are lost.
% Tracking can be 6d or 4d as defined by the input lattice. 
% This is a low level function called by calcDA_raw, which in turn is
% called by the higher level wrapper function "calcDA" or by
% optimzation functions used in MOGA/SOGA
%
%% Usage example
% DAV    = calcDA_border(RING, 0.01, 20, 100, 0.0, 0.0, 0.0, 0.001, 1.01, 0.015, -0.015, 0.007);
%
%% Mandatory input arguments
% r0:     Inital amplitude X to search [m]
% nsteps: Number of radial lines to search     
% etax: dispersion function at tracking position [m]
% res:    Resolution to find the DA
% alpha: factor for dynamic aperture border search.
% nturns: Numbers of turns to track   
% dp: initial momentum or fixed momentum (iof trackiong is 4d
% z0: initial longitudinal position [m]
% xmax: limits of the region where to look for the DA [m]
% xmin: limits of the region where to look for the DA [m]
% ymax: limits of the region where to look for the DA [m]
%
%% Output parameters
% DAV: (nx2) vector containing DA border coordinates [m] 
%                   1 if particle is not lost in nturns
%                   0 if particle is lost in n turns
%

%% History
% Based on model_DA written by M. Munoz
% 2023/09/30 Modified by Pedro F. Tavares
%   allows generic lattice to be given as input
%   allows definition of factor for growing radius (alpha) during dynamic
%   aperture border search
% 2023/10/01: Modified to allow for parallel computation
% 2023/11/19: Modfied  to take dispersion at tracking point as an input (for performance)
% 2023/12/28: Included limits to x and y so as to avoid large aspect
%             ratios/asymmetry
% 2024/03/08: Genrral update and documentation

%% Calculates DA
angle_step=pi/nsteps;
DAV   = zeros(nsteps+1,2);

x0 = etax*dp;
%Check that the chromatic orbit is stable
[~, loss]=ringpass(RING,[x0 0.0 0 0.0 dp z0]',nturns);
if (loss)
    disp('The chromatic closed orbit is not stable. No DA found');
    DAV(1,1)=0;
    DAV(1,2)=0;
else
    parfor i=1:nsteps+1
        look=true; r_stable=0;
        angle = (i-1)*angle_step;
        r=r0;
        while (look)
            x= x0+r*cos(angle);
            y= r*sin(angle);
            if (x>xmax || x<xmin || y>ymax)
                loss=true;
            else
                [~, loss]=ringpass(RING,[x 0.0 y 0.0 dp z0]',nturns);
            end
            
            if (loss)
                if ((r-r_stable) < res)
                    look=false; % we have found the limit of the DA
                    DAV(i,:)=[r_stable*cos(angle),r_stable*sin(angle)] ;
                    r=r_stable;
                else
                    r=(r+r_stable)/2;
                end
            else
                r_stable=r;
                r=r*alpha;
            end
        end
    end
end

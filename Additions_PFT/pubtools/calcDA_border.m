function DA = modelDA_sim_par( RING, r0, nsteps, nturns, dp, z0, etax, res, alpha, xmax, xmin, ymax)
% modelDA( r0, nsteps, nturns, dp, res)
% Evalutes the Dynamic Aperture by tracking
% required arguments
% r0:     Inital amplitude X to search ~ 0.002 m
% nsteps: Number of points to find     ~ 20
% nturns: Numbers of turns to track    ~ 128
% dp:     Energy deviation             ~ 0.0%
% etax: dispersion function at tracking position m
% res:    Resolution to find the DA    ~ 0.0005 m
% alpha: factor for dynamic aperture border search.
% xmax, xmin : positive and negative maximum x
% ymax : maximum y
% Returns the Dynamic aperture border coordinates in DA 


% Written by M. Munoz
% 2023/09/30 Modified by Pedro F. Tavares
%   allows generic lattice to be given as input
%   allows definition of factor for growing radius (alpha) during dynamic
%   aperture border search
%
% 2023/10/01 Modified by Pedro F Tavares to allow for parallel computation
%
% 2023/11/19 Modfied by Pedro F. Tavares to take dispersion at tracking
% point as an input (for performance)
%
% 2023/12/28 Included limits to x and y so as to avoid large aspect
% ratios/asymmetry
%
angle_step=pi/nsteps;
DA   = zeros(nsteps+1,2);
% Evaluate the Chromatic orbit
% twiss=  gettwiss(THERING, 0.0);

% x0=twiss.etax(1)*dp;
% x0=0.0;
x0 = etax*dp;
%Check that the chromatic orbit is stable
[~, loss]=ringpass(RING,[x0 0.0 0 0.0 dp z0]',nturns);
if (loss)
    disp('The chromatic closed orbit is not stable. No DA found');
    DA(1,1)=0;
    DA(1,2)=0;
else
    parfor i=1:nsteps+1
        look=true; r_stable=0;
        angle = (i-1)*angle_step;
        r=r0;
%       fprintf('Tracing step %ld of %ld\n', i, nsteps);
        while (look)
            x= x0+r*cos(angle);
            y= r*sin(angle);
            if (x>xmax || x<xmin || y>ymax)
                loss=true;
            else
                [~, loss]=ringpass(RING,[x 0.0 y 0.0 dp z0]',nturns);
            end
            %fprintf('%s %d %d   \n','Tracked',r, angle);

            if (loss)
                if ((r-r_stable) < res)
                    look=false; % we have found the limit of the DA
                    DA(i,:)=[r_stable*cos(angle),r_stable*sin(angle)] ;
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

function DA = modelDA_sim( RING, r0, nsteps, nturns, dp, res, alpha)
% modelDA( r0, nsteps, nturns, dp, res)
% Evalutes the Dynamic Aperture by tracking
% required arguments
% r0:     Inital amplitude X to search ~ 0.002 m
% nsteps: Number of points to find     ~ 20
% nturns: Numbers of turns to track    ~ 128
% dp:     Energy deviation             ~ 0.0%
% res:    Resolution to find the DA    ~ 0.0005 m
% Returns the Dynamic aperture in DA and the Data for the tracking


% Written by M. Munoz
% 2023/09/30 Modified by Pedro F. Tavares
%   allows generic lattice to be given as input
%   allows definition of factor for growing radius (alpha) during dynamic
%   aperture border search.Other minor changes for performance.

angle_step=pi/nsteps;

angle=0;
r=r0;
DA = zeros(nsteps+1,2);


%Evaluate the Chromatic orbit
% twiss=  gettwiss(THERING, 0.0);

% x0=twiss.etax(1)*dp;
x0=0.0;
%Check that the chromatic orbit is stable
[~, loss]=ringpass(RING,[x0 0.0 0 0.0 dp 0.0]',nturns);
if (loss)
    disp('The chromatic closed orbit is not stable. No DA found');
    DA(1,1)=0;
    DA(1,2)=0;
else
    for i=1:nsteps+1
        look=true; r_stable=0;
%       fprintf('Tracing step %ld of %ld\n', i, nsteps);
        while (look)
            x= x0+r*cos(angle);
            y= r*sin(angle);
            [~, loss]=ringpass(RING,[x 0.0 y 0.0 dp 0.0]',nturns);
            %fprintf('%s %d %d   \n','Tracked',r, angle);

            if (loss)
                if ((r-r_stable) < res)
                    look=false; % we have found the limit of the DA
                    DA(i,1)=r_stable*cos(angle);
                    DA(i,2)=r_stable*sin(angle);
                    r=r_stable;
                else
                    r=(r+r_stable)/2;
                end
            else
                r_stable=r;
                r=r*alpha;
            end
        end
        angle=angle+angle_step;
    end
end


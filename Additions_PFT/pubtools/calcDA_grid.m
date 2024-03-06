function DA = calcDA_grid(RING, X0, Y0, nturns, dp, z0, etax, xmax, xmin, ymax)
%
% Evalutes the Dynamic Aperture by tracking
% required arguments
% RING: the lattice to be evaluated
% X0: a vector of initial horizontal coordinates
% Y0: a vector of initial vertical coordinates
% X0 and Y0 are nx1 column vectors
% initial momenta are assumed to be zero TBC
%
% nturns: Numbers of turns to track    ~ 64
% dp:     Energy deviation             ~ 0.0%
% z0:     LOngitudinal position [m]
% Returns the Dynamic aperture in a nx1 vector :
% 1 if particle is not lost in nturns
% 0 iof particle is lost in n turns
%
%
% 2023/11/06 Written by Pedro F. Tavares
%
%

np   = size(X0,1);
DA   = zeros(np,1);
% Evaluate the Chromatic orbit
% twiss=  gettwiss(THERING, 0.0);

% x0=twiss.etax(1)*dp;
x = etax*dp;
%Check that the chromatic orbit is stable
[~, loss]=ringpass(RING,[x 0.0 0 0.0 dp z0]',nturns);
if (loss)
    disp('The chromatic closed orbit is not stable. No DA found');
    DA =[];
else
    parfor i=1:np
        x = X0(i);
        y = Y0(i);
        if (x>xmax || x< xmin || y >ymax)
            loss=true;
        else
            [~, loss]=ringpass(RING,[x 0.0 y 0.0 dp z0]',nturns);
        end
        DA(i)=not(loss); %
    end
end

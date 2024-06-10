function [Rout, loss, lossinfo] = psTrack(varargin)
% Tracks n particles over a number of turns in a given lattice
%
%% Inputs 
% Mandatory arguments
%
% RING: AT2.0 lattice cell array
% Rin : 6xN matrix: input coordinates of N particles: default=[0 0 0 0 0 0]'
%                   AT 2.0 units . [m rad m rad 1 m]
%
% Optional arguments
%
% S0 : initial azimuthal position to track from; default = 0.0
% S0Tol: tolerance for finding S0 [m]; default = 0.001
% nturns : number of turns
% if any of the two below is nan - no axis limits defined
% Xmax  : horizontal position axis range in plots [mm]; default = nan
% Xmaxp : horizontal angle axis range in plots [mrad]; default  = nan
% if any of the two below is nan - no axis limits defined
% Ymax  : horizontal position axis range in plots [mm]; default = nan
% Ymaxp : horizontal angle axis range in plots [mrad]; default  = nan 
% if any of the two below is nan - no axis limits defined
% deltamax : maximam momentum deviation for plots [%]; default = nan
% zmin     : minimum longitudinal coodinate deviation [mm]; default = nan
% zmax     : maximum longitudinal coodinate deviation [mm]; default = nan
%
% Optional flags
%
% plotx : plots horizontal phase space
% ploty : plots vertical phase space
% plotE : plots vertical phase space
% 6d    : tracking is 6d, otherwise, tracking is 4d. 
% verbose : verbose otput
% 
%% Outputs
% Rout : 6x(N*nturns) matrix with output coordinates of all particles at
%        all turns
% loss : 1xN matrix: 1 if particle lost, 0 if not
% lossinfo: LOSSINFO	1x1 structure with the following fields:
%                lost                 1xN logical vector, indicating lost particles
%                turn                 1xN vector, turn number where the particle is lost
%                element              1xN vector, element number where the particle is lost
%                coordinates_at_loss  6xN array, coordinates at the exit of
%                                     the element where the particle is lost
%                                     (sixth coordinate is inf if particle is lost in a physical aperture)
%                coordinates          6xNxNHIST array, coordinates at the entrance of the
%                                     LHIST elements before the loss
% 
%
%% Usage examples
% Rin = [0.001 0 0 0 0 0]';
% [Rout, loss, lossinfo] = psTrack(RING,Rin, 'S0',528/2,'plotx','nturns',1024,'6d');
% Rout=psTrack(RINGe,Rin,'nturns',500,'verbose','4d','plotE','plotx');
% psTrack(RINGe,Rin,'verbose','6d','plotE','plotx');

%% History
% 2024/03/10 PFT
% 2024/03/10 removed transposing of lattice array before tracking
%            added axis limit control
%            added column to Rout (and plots) with the initial coordinates.
%
%% Input Argument Parsing
[RING, Rin]    = getargs(varargin,[],[0.0 0.0 0.0 0.0 0.0 0.0]');
S0             = getoption(varargin,'S0',0.0);
S0Tol          = getoption(varargin,'S0Tol',0.001);
nturns         = getoption(varargin,'nturns',1024);
Xmax           = getoption(varargin,'Xmax',nan);
Xmaxp          = getoption(varargin,'Xmaxp',nan);
Ymax           = getoption(varargin,'Ymax',nan);
Ymaxp          = getoption(varargin,'Ymaxp',nan);
deltamax       = getoption(varargin,'deltamax',nan);
zmax           = getoption(varargin,'zmax',nan);
zmin           = getoption(varargin,'zmin',nan);
plotxf         = any(strcmpi(varargin,'plotx'));
plotyf         = any(strcmpi(varargin,'ploty'));
plotEf         = any(strcmpi(varargin,'plotE'));
T6df           = any(strcmpi(varargin,'6d'));
verbosef       = any(strcmpi(varargin,'verbose'));

%% Preamble
if (verbosef)
    fprintf('**** \n');
    fprintf('%s Phase Space Tracking %3d  \n', datetime);
end
if (T6df)
    RING = atenable_6d(RING);
else
    RING = atdisable_6d(RING);
end
SPos=findspos(RING,1:length(RING)+1);
Ipos=find(abs(SPos-S0)<=S0Tol);
if isempty(Ipos)
    fprintf('%s Could  not find position to track at %5.3f reset to zero \n', datetime, S0);
    Ipos=1;
end
if (verbosef)
    fprintf('%s Found Initial position to track at %3d \n', datetime, Ipos(1));
end

npart = size(Rin,2);
RING_cy = [RING(Ipos(1):end); RING(1:Ipos(1)-1)];

%% Tracks particles
if (verbosef)
    fprintf('%s Tracking... %3d \n', datetime);
    tic;
end
[Rout,loss,~,lossinfo]=ringpass(RING_cy, Rin, nturns);

Rout = cat(2,Rin,Rout); %%' adds a first
if (verbosef)
    toc;
end

%
%% Plots tracked particles
if(plotxf)
    figure;xlabel('X[mm]');ylabel('Xp[mrad]');
    if (not(isnan(Xmax))&&not(isnan(Xmaxp)))
        xlim([-Xmax Xmax]);ylim([-Xmaxp Xmaxp]);
    end
    grid;hold on;
    for i=1:npart
        plot(Rout(1,i:npart:end)*1000,Rout(2,i:npart:end)*1000,'o'); 
    end
    hold off;

    figure;xlabel('Turn #');ylabel('X[mm]');
    if (not(isnan(Xmax)))
        ylim([-Xmax Xmax]);
    end
    grid;hold on;
    for i=1:npart
        plot(0:nturns,Rout(1,i:npart:end)*1000,'o');
    end
    hold off;
end

if(plotyf)
    figure;xlabel('Y[mm]');ylabel('Yp[mrad]');
    if (not(isnan(Ymax))&&not(isnan(Ymaxp)))
       xlim([-Ymax Ymax]);ylim([-Ymaxp Ymaxp]);
    end
    grid;hold on;
    for i=1:npart
        plot(Rout(3,i:npart:end)*1000,Rout(4,i:npart:end)*1000,'o'); 
    end
    hold off;

    figure;xlabel('Turn #');ylabel('Y[mm]');
    if (not(isnan(Ymax)))
        ylim([-Ymax Ymax]);
    end
    grid;hold on;
    for i=1:npart
        plot(0:nturns,Rout(3,i:npart:end)*1000,'o');
    end
    hold off;
end

if (plotEf)
    figure;xlabel('z[mm]');ylabel('dp[%]');
    if (not(isnan(zmax))&&not(isnan(zmin))&&not(isnan(deltamax)))
       xlim([zmin zmax]);ylim([-deltamax deltamax]);
    end
    grid;hold on;
    for i=1:npart
        plot(Rout(6,i:npart:end)*1000,Rout(5,i:npart:end)*100,'o'); 
    end
    hold off;

    figure;xlabel('Turn #');ylabel('z[mm]');
    if (not(isnan(deltamax)))
        ylim([-deltamax deltamax]);
    end
    grid;hold on;
    for i=1:npart
        plot(0:nturns,Rout(6,i:npart:end)*1000,'o'); 
    end
    hold off;

    figure;xlabel('Turn #');ylabel('dp[%]');
    if (not(isnan(zmax))&&not(isnan(zmin)))
        ylim([zmin zmax]);
    end
    grid;hold on;
    for i=1:npart
        plot(0:nturns,Rout(5,i:npart:end)*100,'o'); 
    end
    hold off;
end


function [Rout, loss] = PSTrack(varargin)
% Tracks n particles over a number of turns in a given lattice
%
%% Mandatory input arguments
% LAT : AT2.0 lattice cell array
% Rin : 6xN matrix: input coordinates of N particles: default=[0 0 0 0 0 0]'
%                   AT 2.0 units . [m rad m rad 1 m]
%
%% Optional arguments
% S0 : initial azimuthmal position to track from; default = 0.0
% S0Tol: tolerance for finding S0 [m]; default = 0.001
% nturns : number of turns
% Xmax  : horizontal position axis range in plots [mm]; default = 7 mm
% Xmaxp : horizontal angle axis range in plots [mrad]; default  = 2.0 mrad
%
%% Optional flags
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
%
% USE:
% [Rout, loss] = PSTrack(RING,[0.001 0 0 0 0 0]', 'S0',528/2,'plotx','nturns',1024,'6d');
%
%% Input Argument Parsing
[LAT, Rin]     = getargs(varargin,[],[0.0 0.0 0.0 0.0 0.0 0.0]');
S0             = getoption(varargin,'S0',0.0);
S0Tol          = getoption(varargin,'S0Tol',0.001);
nturns         = getoption(varargin,'nturns',1024);
Xmax           = getoption(varargin,'Xmax',7.0);
Xmaxp          = getoption(varargin,'Xmaxp',2.0);
plotxf         = any(strcmpi(varargin,'plotx'));
plotyf         = any(strcmpi(varargin,'ploty'));
plotEf         = any(strcmpi(varargin,'plotE'));
T6df           = any(strcmpi(varargin,'6d'));
verbosef       = any(strcmpi(varargin,'verbose'));

%% preamble
if (verbosef)
    fprintf('**** \n');
    fprintf('%s Phase Space Tracking %3d  \n', datetime);
end
if (T6df)
    LAT = atenable_6d(LAT);
else
    LAT = atdisable_6d(LAT);
end
SPos=findspos(LAT,1:length(LAT)+1);
Ipos=find(abs(SPos-S0)<=S0Tol);
if isempty(Ipos)
    fprintf('%s Could  not find position to track at %5.3f reset to zero \n', datetime, S0);
    Ipos=1;
end
if (verbosef)
    fprintf('%s Found Initial position to track at %3d \n', datetime, Ipos(1));
end

npart = size(Rin,2);
LAT_cy = [LAT(Ipos(1):end); LAT(1:Ipos(1)-1)]';
if (verbosef)
    fprintf('%s Tracking... %3d \n', datetime);
    tic;
end
[Rout,loss]=ringpass(LAT_cy, Rin, nturns);

if (verbosef)
    toc;
end

if(plotxf)
    figure;xlabel('X[mm]');ylabel('Xp[mrad]');xlim([-Xmax Xmax]);ylim([-Xmaxp Xmaxp]);grid;hold;
    for i=1:npart
        plot(Rout(1,i:npart:end)*1000,Rout(2,i:npart:end)*1000,'o'); 
    end
    figure;xlabel('Turn #');ylabel('X[mm]');ylim([-Xmax Xmax]);grid;hold;
    for i=1:npart
        plot(Rout(1,i:npart:end)*1000,'o');
    end
end

if(plotyf)
    figure;xlabel('Y[mm]');ylabel('Yp[mrad]');hold;
    for i=1:npart
        plot(Rout(3,i:npart:end)*1000,Rout(4,i:npart:end)*1000,'o'); 
    end
    
end

if (plotEf)
    figure;xlabel('z[mm]');ylabel('dp[%]');hold;
    for i=1:npart
        plot(Rout(6,i:npart:end)*1000,Rout(5,i:npart:end)*100,'o'); 
    end
    figure;xlabel('Turn #');ylabel('z[mm]');hold;
    for i=1:npart
        plot(Rout(6,i:npart:end)*1000,'o'); 
    end
    figure;xlabel('Turn #');ylabel('dp[%]');hold;
    for i=1:npart
        plot(Rout(5,i:npart:end)*100,'o'); 
    end
end


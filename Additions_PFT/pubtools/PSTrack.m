function [Rout, loss] = PSTrack(varargin)
% Tracks n particles over a nuebr of ituns in a given lattice
%
%% Inputs
% Mandatory input arguments
% LAT : AT2.0 lattice cell array
%
% Optional arguments
% Rin : 6xN matrix: input coordinates of N particles - default = zero.
% nturns : number of turns
% 'plot' if present, plots phase space
%
%% Input Argument Parsing
LAT            = getargs(varargin,[]);
Rin            = getoption(varargin,'Rin',[0.001 0.0 0.0 0.0 0.0 0.0]');
nturns         = getoption(varargin,'nturns',1024);
Xmax           = getoption(varargin,'Xmax',7.0);
Xmaxp          = getoption(varargin,'Xmaxp',2.0);
plotxf         = any(strcmpi(varargin,'plotx'));
plotyf         = any(strcmpi(varargin,'ploty'));
plotEf         = any(strcmpi(varargin,'plotE'));
T6df           = any(strcmpi(varargin,'6d'));



%% preamble
if (T6df)
    LAT = atenable_6d(LAT);
else
    LAT = atdisable_6d(LAT);
end
atsummary(LAT)
npart = size(Rin,2);
[Rout,loss]=ringpass(LAT, Rin, nturns);

if(plotxf)
    figure;xlabel('X[mm]');ylabel('Xp[mrad]');hold;
    for i=1:npart
        plot(Rout(1,i:npart:end)*1000,Rout(2,i:npart:end)*1000,'o'); xlim([-Xmax Xmax]); ylim([-Xmaxp Xmaxp]);grid;
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


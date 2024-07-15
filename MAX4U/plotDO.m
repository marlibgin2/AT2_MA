function plotDO(LattStruct,varargin)
% plots the design Orbit in lab coordinates and the deviation wrt to
% a reference orbit.
% plots the trajectory along the centre of the magnets. This may differ
% from the above in the case of reverse bends implemented as offset
% quadrupoles
%   
%% Inputs
% Mandatory argument
% LattStru : structure generated with the cLatt function cntainign the 
%            design orbit for a lattuice and a reference lattice
%
% Optional flafs
% 'do' : plots design orbit and reference orbit, if available
% 'dev': plots deviation between design orbit and reference orbit
%
%
%
%% Usage examples
% plotDO(m4_standard);

%% History
% PFT 2024/06/13, first version
% PFT 2024/06/30 : added latice title
% PFT 2024/07/09 : added magnet plots of magnet centres
%% Input argument parsing
%
plotdof      = any(strcmpi(varargin,'do'));
plotdevf     = any(strcmpi(varargin,'dev'));
plotmagf     = any(strcmpi(varargin,'mag'));
plotmagdevf  = any(strcmpi(varargin,'magdev'));


%% Plots DO
lattname= LattStruct.Lattice_Name;
x2d     = LattStruct.LattData.DesignOrbit.x2d;
y2d     = LattStruct.LattData.DesignOrbit.y2d;
x2d_ref = LattStruct.LattData.DesignOrbit.x2d_ref;
y2d_ref = LattStruct.LattData.DesignOrbit.y2d_ref;
s2d = LattStruct.LattData.DesignOrbit.s2d;
dev = LattStruct.LattData.DesignOrbit.Deviation;

if (plotdof)
    figure;plot(x2d, y2d, '-ob');hold; plot(x2d_ref, y2d_ref, '-or');
    xlabel('X[m]');ylabel('Y[m]');legend({'design orbit';'reference orbit'});
    grid on; 
    title(lattname);
end

if(plotdevf)
    figure;plot(x2d, dev*1000, '-o'); 
    xlabel('X[m]');ylabel('dZ[mm]');
    grid on;
    title(strcat(lattname,' Design Orbit Deviation'));
end


%% Plots magnet centres
x2d_mag  = LattStruct.LattData.MagCentres.x2d;
y2d_mag  = LattStruct.LattData.MagCentres.y2d;
magdev = LattStruct.LattData.MagCentres.Deviation;

if(plotmagf)
    figure;plot(x2d_mag, y2d_mag, '-ob');hold; plot(x2d, y2d, '-or');
    xlabel('X[m]');ylabel('Y[m]');legend({'Magnet Centres';'Design Orbit'});
    grid on; 
    title(lattname);
end

if(plotmagdevf)
    figure;plot(x2d, magdev*1000, '-o'); 
    xlabel('X[m]');ylabel('dZ[mm]');
    grid on;
    title(strcat(lattname, ' Deviation Magnet Centres'));
end


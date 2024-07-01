function plotDO(LattStruct,varargin)
% plots the Deisgn Orbit in lab coordinates and the deviation wrt to
% a reference orbit
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
%% Input argument parsing
%
plotdof  = any(strcmpi(varargin,'do'));
plotdevf = any(strcmpi(varargin,'dev'));


%% Plots DO
lattname=LattStruct.Lattice_Name;
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
    title(lattname);
end




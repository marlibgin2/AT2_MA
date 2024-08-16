function plotGO(LattStruct,varargin)
% plots the design Orbit in lab coordinates and the deviation wrt to
% a reference orbit.
% plots the trajectory along the centre of the magnets. This may differ
% from the above in the case of reverse bends implemented as offset
% quadrupoles
%   
%% Inputs
% Mandatory argument
% LattStruct : structure generated with the cLatt function containing the 
%            design orbit for a lattice and a reference lattice
%
% Optional flags
% 'do'          : plots design orbit and reference orbit, if available
% 'dodev'       : plots deviation between design orbit and reference orbit
% 'showBPMs'    : indicates BPM positions in orbit deviation plot
% 'magce'       : plots magnet centres
% 'magcedev '   : plots deviation of magnet centres to reference
% 'magap'       : plots magnet apertures
% 'maggeo'      : plots limits of magnet aperture
% 'maggeodev'   : plots deviations of limits of magnet aperture wrt to
%                 reference
% 'magceTobeam' : plots magnet centre-to-beam deviation
% 'chace'       : plots chamber centres
% 'chacedev'    : plots deviation of chamber centres to reference
% 'chaap'       : plots chamber apertures
% 'chageo'      : plots limits of chamber aperture
% 'chageodev'   : plots deviations of limits of chamber aperture wrt to
% 'chaceTobeam' : plots vertical distance from  chamber centre to the beam.
%               reference 
% 'effap'       : plots effective aperture

% 'all'         : do all plots
% 'dos'         : do all design orbit plots
% 'mags'        : do all magnet plots
% 'chas'        : do all chamber plots
%
% 'env'         : plots magnet walls, chamber walls and beam
%% Usage examples
% plotGO(m4_standard);

%% History
% PFT 2024/06/13, first version
% PFT 2024/06/30 : added latice title
% PFT 2024/07/09 : added magnet plots of magnet centres
% PFT 2024/07/18 : added plotting of magnet/chamber geometry
% PFT 2024/07/22 : changed name from plotDO to plotGO 
% PFT 2024/08/16 : added option to indicate BPM positions in the orbit
%                  deviation plot
%% Input argument parsing
%
plotdof         = any(strcmpi(varargin,'do'));
plotdodevf      = any(strcmpi(varargin,'dodev'));
showBPMsf       = any(strcmpi(varargin,'showBPMs'));

plotmagcef      = any(strcmpi(varargin,'magce'));
plotmagcedevf   = any(strcmpi(varargin,'magcedev'));
plotmagapf      = any(strcmpi(varargin,'magap'));
plotmaggeof     = any(strcmpi(varargin,'maggeo'));
plotmaggeodevf  = any(strcmpi(varargin,'maggeodev'));
plotmagceTobeamf= any(strcmpi(varargin,'magceTobeam'));

plotchacef      = any(strcmpi(varargin,'chace'));
plotchacedevf   = any(strcmpi(varargin,'chacedev'));
plotchaapf      = any(strcmpi(varargin,'chaap'));
plotchageof     = any(strcmpi(varargin,'chageo'));
plotchageodevf  = any(strcmpi(varargin,'chageodev'));
plotchaceTobeamf= any(strcmpi(varargin,'chaceTobeam'));
plotefffapf     = any(strcmpi(varargin,'effap'));

plotallf        = any(strcmpi(varargin,'all'));
plotdosf        = any(strcmpi(varargin,'dos'));
plotmagsf       = any(strcmpi(varargin,'mags'));
plotchasf       = any(strcmpi(varargin,'chas'));

plotenvf        = any(strcmpi(varargin,'env'));


%% preamble
lattname      = LattStruct.Lattice_Name;
geometry      = LattStruct.LattData.geometry;
ACHRO         = LattStruct.ACHROMAT;
x2d           = geometry.DesignOrbit.x2d;
y2d           = geometry.DesignOrbit.y2d;
s2d           = geometry.DesignOrbit.s2d;
x2d_ref       = geometry.ref.DesignOrbit.x2d;
y2d_ref       = geometry.ref.DesignOrbit.y2d;
orbdev        = geometry.DesignOrbit.Deviation;


magHAperture  = geometry.Magnets.HAperture;
x2d_magce     = geometry.Magnets.Centre.x2d;
y2d_magce     = geometry.Magnets.Centre.y2d;
magdevce      = geometry.Magnets.Centre.Deviation;

magceTobeam   = geometry.Magnets.magceTobeam;

x2d_magup     = geometry.Magnets.Walls.x2d_up;
y2d_magup     = geometry.Magnets.Walls.y2d_up;
magdevup      = geometry.Magnets.Walls.dev_up;

x2d_magdown  = geometry.Magnets.Walls.x2d_down;
y2d_magdown  = geometry.Magnets.Walls.y2d_down;
magdevdown   = geometry.Magnets.Walls.dev_down;


chaHAperture = geometry.Chambers.HAperture;
x2d_chace    = geometry.Chambers.Centre.x2d;
y2d_chace    = geometry.Chambers.Centre.y2d;
chadevce     = geometry.Chambers.Centre.Deviation;

chaceTobeam  = geometry.Chambers.chaceTobeam;
effectiveAperture = geometry.Chambers.effectiveAperture;

x2d_chaup    = geometry.Chambers.Walls.x2d_up;
y2d_chaup    = geometry.Chambers.Walls.y2d_up;
chadevup     = geometry.Chambers.Walls.dev_up;

x2d_chadown  = geometry.Chambers.Walls.x2d_down;
y2d_chadown  = geometry.Chambers.Walls.y2d_down;
chadevdown   = geometry.Chambers.Walls.dev_down;


iBPM = findcells(ACHRO,'FamName','BPM');
if (isempty(iBPM))
    iBPM=findcells(ACHRO,'FamName','mon');
end

if (isempty(iBPM))
    fprintf('%s Error in calcOrb; no BPMs in Lattice to plot orbit, aborting... \n', datetime);
end
sBPM = findspos(ACHRO,iBPM);
nBPMs=numel(sBPM);
if (nBPMs>0)
    [~, ia, ~] = unique(s2d);
    x2d_u = x2d(ia);
    s2d_u = s2d(ia);
    orbdev_u = orbdev(ia);
    xBPM=zeros(nBPMs,1);
    ydevBPM=zeros(nBPMs,1);
    for i=1:nBPMs
        xBPM(i)=interp1(s2d_u,x2d_u,sBPM(i));
        ydevBPM(i)=interp1(x2d_u,orbdev_u*1000,xBPM(i));
    end
end

%% Plots design orbit 
if (plotdof||plotdosf||plotallf)
    figure;plot(x2d, y2d, '-ob');hold on; plot(x2d_ref, y2d_ref, '-or');
    xlabel('X[m]');ylabel('Y[m]');legend({'design orbit';'reference orbit'});
    grid on; 
    title(lattname);
end

%% Plots design orbit deviation from reference
if(plotdodevf||plotdosf||plotallf)
    figure;plot(x2d, orbdev*1000, '-o');hold on;
    if ((nBPMs>0)&&(showBPMsf))
        plot(xBPM,ydevBPM, 'rs', MarkerSize=16);
    end

    xlabel('X[m]');ylabel('dZ[mm]');
    grid on; 
    title(strcat(lattname,' Design Orbit Deviation'));
    if ((nBPMs>0)&&(showBPMsf))
        legend({'orbit deviation';'BPMs'});
    end
end


%% Plots magnet centres
if(plotmagcef||plotmagsf||plotallf)
    figure;plot(x2d_magce, y2d_magce, '-ob');hold on; plot(x2d, y2d, '-or');
    xlabel('X[m]');ylabel('Y[m]');legend({'Magnet Centres';'Design Orbit'});
    grid on; 
    title(lattname);
end

%% Plots magnet centre deviations from reference
if(plotmagcedevf||plotmagsf||plotallf)
    figure;plot(x2d, magdevce*1000, '-o'); 
    xlabel('X[m]');ylabel('dZ[mm]');
    grid on;
    title(strcat(lattname, ' Deviation Magnet Centres'));
end
%% Plots magnet centre to beam deviation
if(plotmagceTobeamf||plotmagsf||plotallf)
    figure;plot(x2d, magceTobeam*1000, '-b'); 
    xlabel('X[m]');ylabel('dZ[mm]');
    grid on;
    title(strcat(lattname, ' Magnet-centre-to-beam deviation'));
end
%% Plots magnet apertures

if(plotmagapf||plotmagsf||plotallf)
    figure;plot(x2d, magHAperture*1000, '-ob');
    xlabel('X[m]');ylabel('Horizontal Magnet Aperture[mm]');
    grid on; 
    title(lattname);
end


%% Plots magnet geometry

if(plotmaggeof||plotmagsf||plotallf)
    figure;plot(x2d_magup, y2d_magup, '-ob'); hold on;
    plot(x2d_magdown, y2d_magdown, '-or');
    plot(x2d_magce, y2d_magce, '--ok');
    plot(x2d,y2d,'-g')

    xlabel('X[m]');ylabel('Y[m]');legend({'Outer';'Inner';'Centre';'beam'});
    grid on; 
    title(strcat(lattname, {' Magnet Geometry'}));
end


%% Plots magnet geometry deviation from reference
if(plotmaggeodevf||plotmagsf||plotallf)
    figure;plot(x2d, magdevup*1000, '-b'); hold on;
    plot(x2d, magdevdown*1000, 'or');
    xlabel('X[m]');ylabel('Deviation [mm]');legend({'Outer';'Inner'});
    grid on; 
    title(strcat(lattname,{' Magnet Aperture'}));
end
%% Plots chamber centres
if(plotchacef||plotchasf||plotallf)
    figure;plot(x2d_chace, y2d_chace, '-ob');hold on; plot(x2d, y2d, '-r');
    xlabel('X[m]');ylabel('Y[m]');legend({'Chamber Centre';'Design Orbit'});
    grid on; 
    title(lattname);
end

%% Plots chamber centre deviations from reference
if(plotchacedevf||plotchasf||plotallf)
    figure;
    if (numel(x2d)==numel(chadevce))
        plot(x2d, chadevce*1000, '-o'); 
    else
        plot(x2d_ref, chadevce*1000, '-o'); 
    end


    xlabel('X[m]');ylabel('dZ[mm]');
    grid on;
    title(strcat(lattname, ' Chamber Centre Deviation'));
end
%% Plots chamber centre to beam deviation
if(plotchaceTobeamf||plotchasf||plotallf)
    figure;
    if (numel(x2d)==numel(chaceTobeam))
        plot(x2d, chaceTobeam*1000, '-o'); 
    else
        plot(x2d_ref, chaceTobeam*1000, '-o'); 
    end

    xlabel('X[m]');ylabel('dZ[mm]');
    grid on;
    title(strcat(lattname, ' Chamber-centre-to-beam deviation'));
end
%% Plots chamber aperture
if (plotchaapf||plotchasf||plotallf)
    figure;
    if (numel(x2d)==numel(chaHAperture))
        plot(x2d, chaHAperture*1000, '-b');
    else
        plot(x2d_ref, chaHAperture*1000, '-b'); 
    end
    hold on;
    xlabel('X[m]');ylabel('Horizontal Aperture[mm]');
    grid on; 
    title(lattname);
end
%% Plots Effective aperture
if (plotefffapf||plotchasf||plotallf)
    figure;
    if (numel(x2d)==numel(effectiveAperture))
        plot(x2d, effectiveAperture*1000, '-b');
    else
        plot(x2d_ref, effectiveAperture*1000, '-b'); 
    end
    hold on;
    xlabel('X[m]');ylabel('Horizontal Chamber Aperture[mm]');
    grid on; 
    title(strcat(lattname, {' Effective Aperture'}));
end

%% Plots chamber geometry
if (plotchageof||plotchasf||plotallf)
    figure;plot(x2d_chaup, y2d_chaup, '-ob'); hold on;
    plot(x2d_chadown, y2d_chadown, '-or');
    plot(x2d_chace,y2d_chace,'--k');
    plot(x2d,y2d,'g');

    xlabel('X[m]');ylabel('Y[m]');legend({'Outer';'Inner';'Centre';'beam'});
    grid on; 
    title(strcat(lattname, {' Chamber Geometry'}));
end


%% Plots chamber geometry deviation from reference
if(plotchageodevf||plotchasf||plotallf)
    figure;
    if (numel(x2d)==numel(chadevup))
        x2dl=x2d;
    else
        x2dl=x2d_ref;
    end
    plot(x2dl, chadevup*1000, '-b'); hold on;
    plot(x2dl, chadevdown*1000, 'or');
    xlabel('X[m]');ylabel('Deviation [mm]');legend({'Outer';'Inner'});
    grid on; 
    title(strcat(lattname,{' Chamber Geometry Deviation'}));
end
%% Plots full envelope, i.e. magnet walls, chamber walls and beam
if(plotenvf ||plotchasf||plotmagsf||plotallf)
    figure;plot(x2d_chaup, y2d_chaup, '--b'); hold on;
    plot(x2d_chadown, y2d_chadown, '--r');
    
    plot(x2d_magup, y2d_magup, '-b'); hold on;
    plot(x2d_magdown, y2d_magdown, '-r');

    plot(x2d,y2d,'g');

    xlabel('X[m]');ylabel('Y[m]');legend({'Outer Chamber';'Inner Chamber';'Outer Magnet';'Inner Magnet';'beam'});
    grid on; 
    title(strcat(lattname, {' Full Envelope'}));
end

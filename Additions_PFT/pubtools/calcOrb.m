function [orb0,varargout] = calcOrb(varargin)
% Calculates, plots and corrects the closed orbit
% 
% This is a higher level wrapper function
% 
%% Usage examples
% [orb0, orb, RINGc] = calcOrb(RING,'plot');
%
%% Mandatory input arguments
% RING : AT2 lattice array
%
%% Optional input parameters
% None
%
%% Optional flags
% correct: corrects the orbit
% plot : plots the orbit
% verbose: produces verbose output
%
%% Output parameters
% orb0: (nx2) array: (X,Y) Initial orbit  [m] DAoptions: Structure with options used in the calculation
% RINGc: corrected ring (only if correction is asked for)
% orb: (nx2) array: (X,Y) Corrected orbit [m] (only if correction is asked
% 

%% History
% PFT 2024/03/02
%
%% input argument parsing
RING           = getargs(varargin,[]);
plotf          = any(strcmpi(varargin,'plot'));
correctf       = any(strcmpi(varargin,'correct'));
verbosef       = any(strcmpi(varargin,'verbose'));

%% Calculates the closed orbit
% View the orbit, including BPM errors
setoption('WarningDp6D',false); % avoids warning messages
iBPM = findcells(RING,'FamName','BPM');
if (isempty(iBPM))
    fprintf('%s Error in calcOrb; no BPMs in Lattice to plot orbit, aborting... \n', datetime);
    orb0=nan;
    varargout{1}=nan;
    varargout{2}=nan;
    return
end
sBPM = findspos(RING,iBPM);
orb0 = findorbit6Err(RING,iBPM);

if (plotf)
    figure; plot(sBPM,1e3*orb0([1 3],:)); xlim([0 528]);
            xlabel('s [m]'); ylabel('x,y [mm]');grid;legend('X','Y');
    hold on;
end

% Correct the orbit
if (correctf)
%    RINGc = atcorrectorbit(RING,[],[],[],[],[140 120; 160 140; 180 160; ...
%                           ones(10,1)*[200 180]],[true true],0.75,...
%                           [],[],[0.38, 0.38]*1e-3,verbosef);
    RINGc = atcorrectorbit(RING,[],[],[],[],[],[true true],0.75,...
                           [],[],[0.38, 0.38]*1e-3,verbosef);
    varargout{1}=RINGc;
    % Calculate the new orbit and plot in the former figure
    orb = findorbit6Err(RINGc,iBPM);
    varargout{2}=orb;
    if (plotf)
        figure;plot(sBPM,1e6*orb([1 3],:));xlim([0 528]); 
               xlabel('s [m]'); ylabel('x,y [µm]');grid;legend('X','Y');
    end
else
    varargout{1}=nan;
    varargout{2}=nan;
end






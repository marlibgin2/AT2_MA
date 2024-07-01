function [RINGc,orb0,orb] = calcOrb(varargin)
% Calculates, plots and corrects the closed orbit
% 
% This is a higher level wrapper function
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
% RINGc: corrected ring (if correction is not asked for, this is the same as the input lattice)
% orb0: (nx2) array: (X,Y) Initial orbit  [m] 
% orb: (nx2) array: (X,Y) Corrected orbit [m] (only if correction is done)
% 
%% Usage examples
% [RINGc, orb0, orb] = calcOrb(RING,'plot','correct');
% [~, orb0, orb]     = calcOrb(RING,'plot');
% RINGc = calcOrb(RING,'plot','correct');
% calcOrb(RING,'plot');

%% History
% PFT 2024/03/02, first version
% PFT 2024/03/10: Bug fix to get the correct output lattice
% PFT 2024/06/30: added search for corrector indices based on family names
%
%% Input argument parsing
RING           = getargs(varargin,[]);
plotf          = any(strcmpi(varargin,'plot'));
correctf       = any(strcmpi(varargin,'correct'));
verbosef       = any(strcmpi(varargin,'verbose'));

%% Calculates the closed orbit
% View the orbit, including BPM errors
setoption('WarningDp6D',false); % avoids warning messages
iBPM = findcells(RING,'FamName','BPM');
RINGc = RING;
if (isempty(iBPM))
    fprintf('%s Error in calcOrb; no BPMs in Lattice to plot orbit, aborting... \n', datetime);
    orb0=nan;
    orb=nan;
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
    indHCor=find(atgetcells(RING,'iscorH','H'));
    if (isempty(indHCor))
        indHCor=findcells(RING,'FamName','ch');
    end
    indVCor=find(atgetcells(RING,'iscorV','V'));
    if (isempty(indHCor))
        indVCor=findcells(RING,'FamName','ch');
    end

%    RINGc = atcorrectorbit(RING,[],[],[],[],[140 120; 160 140; 180 160; ...
%                           ones(10,1)*[200 180]],[true true],0.75,...
%                           [],[],[0.38, 0.38]*1e-3,verbosef);
     RINGc = atcorrectorbit(RING,iBPM,indHCor,indVCor,[],[],[true true],0.75,...
                           [],[],[0.38, 0.38]*1e-3,verbosef);
    % Calculate the new orbit and plot in the former figure
    orb = findorbit6Err(RINGc,iBPM);
    if (plotf)
        figure;plot(sBPM,1e6*orb([1 3],:));xlim([0 528]); 
               xlabel('s [m]'); ylabel('x,y [µm]');grid;legend('X','Y');
    end   
end






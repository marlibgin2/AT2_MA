function [orb0, orb, RINGc] = calcOrb(varargin)
% Calculates, Plots and corrects the closed orbit
% 
RING           = getargs(varargin,[]);
plotf          = any(strcmpi(varargin,'plot'));
correctf       = any(strcmpi(varargin,'correct'));
verbosef       = any(strcmpi(varargin,'verbose'));

% View the orbit, including BPM errors
iBPM = findcells(RING,'FamName','BPM');
sBPM = findspos(RING,iBPM);
orb0 = findorbit6Err(RING,iBPM);

if (plotf)
    figure; plot(sBPM,1e6*orb0([1 3],:)); xlim([0 528]); xlabel('s [m]'); ylabel('x,y [µm]')
    hold on;
end

% Correct the orbit
if (correctf)
    RINGc = atcorrectorbit(RING,[],[],[],[],[140 120; 160 140; 180 160; ones(10,1)*[200 180]],[true true],0.75,[],[],[0.38, 0.38]*1e-3,verbosef);
    % Calculate the new orbit and plot in the former figure
    orb = findorbit6Err(RINGc,iBPM);
    if (plotf)
        plot(sBPM,1e6*orb([1 3],:));
    end
end





end
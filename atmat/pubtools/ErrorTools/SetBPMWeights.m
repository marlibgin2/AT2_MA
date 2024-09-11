function ringW=SetBPMWeights(ring)
% In order to avoid issues with orbit correction routines, which all seem
% to assume that we have as many vertical correctors as we do BPMs, we need
% to add weight information to the BPMs (and make sure correction routines
% make use of it).
%
    ringW=ring;
    I = findcells(ringW,'FamName','BPM');
    if (isempty(I))
        I=findcells(ringW,'FamName','mon');
    end
    
    if (not(isempty(I)))
        ringW{I(1)}.Weight = [0 0]; I = I(2:end);
        for n = 1:numel(I)
            if mod(n-2,10) == 0
                ringW{I(n)}.Weight = [1 1e-3];
            else
                ringW{I(n)}.Weight = [1 1];
            end
        end
    else
        fprintf('%s SetBPMWeigths: Warning - no BPMs found. \n', datetime);
    end
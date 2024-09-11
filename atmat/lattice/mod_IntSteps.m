function RING = mod_IntSteps(RING,nIStep_DP,nIStep_QP,nIStep_SP,nIStep_OP)
%% Check all the elements for the Number of integration steps defined in lattice file 
%  It sets numIntSteps to the desired input value for dipole, quadrupole,
%  sextupole and octupole. It will set numIntSteps= 1 for all other elements 
%% Input
%  RING: AT2 lattice array 
%  nIStep_DP: Number of integrtion steps for dipoles
%  nIStep_QP: Number of integrtion steps for Quadrupoles
%  nIStep_SP: Number of integrtion steps for Sextupoles
%  nIStep_OP: Number of integrtion steps for Octupoles
%% Usage examples
% RING = mod_IntSteps(RING,nIStep_DP,nIStep_QP,nIStep_SP,nIStep_OP)
% RING = mod_IntSteps(RING,10,3,2,1);

%% History
% S. Jena  2024/07/16

% Loop through all elements in the lattice
for i = 1:length(RING)
    element = RING{i};
    
    numIntSteps = 1;  % Default for elements not specified

    if isfield(element, 'PassMethod')
        % Check if the element is a dipole
        if strcmp(element.PassMethod, 'BndMPoleSymplectic4Pass')
            numIntSteps = nIStep_DP; %10
        elseif strcmp(element.PassMethod, 'StrMPoleSymplectic4Pass') && isfield(element, 'PolynomB')
            % Check the PolynomB field to distinguish between quadrupoles, sextupoles, and octupoles
            if length(element.PolynomB) > 1 && element.PolynomB(2) ~= 0
                numIntSteps = nIStep_QP;  % Quadrupole -3
            elseif length(element.PolynomB) > 2 && element.PolynomB(3) ~= 0
                numIntSteps = nIStep_SP;  % Sextupole -2
            elseif length(element.PolynomB) > 3 && element.PolynomB(4) ~= 0
                numIntSteps = nIStep_OP;  % Octupole -1
            end
        end
    end

    % Set the NumIntSteps field for the element
    element.NumIntSteps = numIntSteps;
    RING{i} = element;
end

% Save the modified lattice to a new file
% save('RINGn.mat', 'RING');

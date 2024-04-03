function EM = errormodel_example
% ERRORMODEL_EXAMPLE contains an example (unrealistic) demonstration error model 
%
% NOTES
% 1. Intended usage is with applyErrorModel, which will deploy the defined
%    model on a lattice.
%
% 2. The return structure has a number of fields corresponding to the
%    'units' on which the errors should be deployed. E.g. girders, or
%    magnets (which may be sliced in the lattice).
%
% 3. For magnets, the errors can be defined based on AT2 class (e.g. 'Bend'
%    or 'Quadrupole') or FamilyName. Specific trumps general, so if errors
%    have been defined for a quadrupole named 'QF' -only- those errors will
%    be applied. Generic quadrupole errors will not be applied in that
%    case. This allows the use of generic quad errors that can then be
%    selectively more accurate for certain magnets. 
%      NB #1! If errors are defined for a specific FamName, even if only a
%       single misalignment, ALL applicable class errors will be ignored.
%       So if you do define e.g. a field error for a quad family 'QF'
%       you must also define misalignment errors for the 'QF' family
%       explicitly! 
%      NB #2! It is perfectly possible to define errors to define to
%       multiple FamNames. Make sure the ID field has a cell array
%       containing all FamNames to which the error is applicable (e.g.
%       {'R1','R2','R3'}). This is particularly useful to distinguish
%       between ordinary bends and reverse bends.
%    Note that non-magnet classes may be used as well, regardless of
%    whether it makes sense or not. 
%
% 4. For girders, it is recommended to always have a 'Baseline' alignment
%    accuracy defined. Note that it is possible to define different errors
%    depending on the girder type. To do this, ensure that the girder
%    marker elements GS and GE in the lattice have an additional field
%    named 'GirderType'. It is then possible to specify different errors
%    for different girders, e.g. 'UnitCell5', in the error model.
%      NB #1! As for the magnets, specific trumps generic. If the
%      GirderType in the lattice matches the error model ID, only the
%      errors defined for that girder will be deployed and 'Baseline' errors
%      will be ignored.
%
% See also applyErrorModel, getMagGroupsFromGirderIndex, getMagGroupsFromMagNum


%% Girder error definitions
% Error model used is first based on girder name, which will be a


% Generic girder errors
Girder{1}.ID = 'Baseline';

Girder{1}.Systematic{1} = struct( ...
    'Pitch', 0, ...
    'Yaw', 0, ...
    'Roll', 0, ...
    'Sway', 0, ...
    'Heave', 0, ...
    'Surge', 0);
Girder{1}.Random{1} =     struct( ...
    'Pitch', 0e-6, ...
    'Yaw', 0e-6, ...
    'Roll', 0e-3, ...
    'Sway', 0e-6, ...
    'Heave', 0e-6, ...
    'Surge', 0);


% Generic girder errors
Girder{1}.ID = 'U1';        % This is matched against the girder marker 'GirderType' field. See help text. 

Girder{1}.Systematic{1} = struct( ...
    'Pitch', 0, ...
    'Yaw', 0, ...
    'Roll', 0, ...
    'Sway', 0, ...
    'Heave', 100e-6, ...
    'Surge', 0);
Girder{1}.Random{1} =     struct( ...
    'Pitch', 0e-6, ...
    'Yaw', 0e-6, ...
    'Roll', 0.1e-3, ...
    'Sway', 0e-6, ...
    'Heave', 0e-6, ...
    'Surge', 0);



%% MAGNET ERROR DEFINITIONS

% Exact same magnet errors are applied to all slices with the same MagNum
% Error magnitude deployed to each slice are scaled with the base strength
% of the slice (NB, not true for e.g. dipoles which get sextupole fringe
% fields. Use systematics and PolynomB of Nx4 size where N=nbr of slices)

% BASELINE MAGNET ALIGNMENT
DDR_ChallengingMachiningTolerances.Systematic = struct( ...
    'Pitch', 0, ...
    'Yaw', 0, ...
    'Roll', 0, ...
    'Sway', 0, ...
    'Heave', 0, ...
    'Surge', 0);
DDR_ChallengingMachiningTolerances.Random = struct( ...
    'Pitch', 0, ...
    'Yaw', 0, ...
    'Roll', 0.1e-3, ...
    'Sway', 0e-6, ...
    'Heave', 0e-6, ...
    'Surge', 0);
DDR_ChallingingBPMCalibrationAccuracy.Random = struct(...
    'Roll', 0e-3, ...
    'Sway', 0e-6, ...
    'Heave', 0e-6);


% GENERIC BEND MODEL
Magnet{1}.ID = 'Bend';      % May be either atclass as here, or FamName. Latter takes precedence if available.

Magnet{end}.Systematic{1}     =   DDR_ChallengingMachiningTolerances.Systematic;
Magnet{end}.Random{1}         =   DDR_ChallengingMachiningTolerances.Random;

% Random errors, incl. multipoles
PolB(1,[2 3 4 6 10]) = [2.5, 2.8, 1.9, 1.3, 0.3]*1e-4;
PolA(1,[3,4]) = [2.9, 1.4]*1e-4;
Magnet{end}.Random{2} = struct( ...
    'PolynomB', PolB, ...     % All values are RMS, and specify variation of and relative to the main component
    'PolynomA', PolA, ...     % All values are RMS, and specify variation of and relative to the main component
    'Scaling', 1.001);        % Multiplicative factor, useful for current errors

% Systematic errors, incl. multipoles
PolB(1,[6 10 14]) = [0.5, 0.5, 0.1]*1e-4;
PolA = [0,0,0,0];
Magnet{end}.Systematic{2} = struct( ...
    'PolynomB', PolB, ...     % All values are RMS
    'PolynomA', PolA, ...     % All values are RMS
    'Scaling', 1);            % Multiplicative factor, useful for current errors



% GENERIC QUADRUPOLE MODEL
Magnet{end+1}.ID = 'Quadrupole';    

Magnet{end}.Systematic{1}     =   DDR_ChallengingMachiningTolerances.Systematic;
Magnet{end}.Random{1}         =   DDR_ChallengingMachiningTolerances.Random;

% Random errors, incl. multipoles
PolB(1,[2 3 4 6 10]) = [2.5, 2.8, 1.9, 1.3, 0.3]*1e-4;
PolA(1,[3,4]) = [2.9, 1.4]*1e-4;
Magnet{end}.Random{2} = struct( ...
    'PolynomB', PolB, ...     % All values are RMS, and specify variation of and relative to the main component
    'PolynomA', PolA, ...     % All values are RMS, and specify variation of and relative to the main component
    'Scaling', 1.001);            % Multiplicative factor, useful for current errors

% Systematic errors, incl. multipoles
PolB(1,[6 10 14]) = [0.5, 0.5, 0.1]*1e-4;
PolA = [0,0,0,0];
Magnet{end}.Systematic{2} = struct( ...
    'PolynomB', PolB, ...     % All values are RMS
    'PolynomA', PolA, ...     % All values are RMS
    'Scaling', 1);            % Multiplicative factor, useful for current errors



% GENERIC SEXTUPOLE MODEL
Magnet{end+1}.ID = 'Sextupole';    % May be either atclass or FamName. Latter takes precedence if available.

Magnet{end}.Systematic{1}   =   DDR_ChallengingMachiningTolerances.Systematic;
Magnet{end}.Random{1}       =   DDR_ChallengingMachiningTolerances.Random;

PolB = zeros(1,21); PolB([3 4 5 9 15]) = [5.0, 5.2, 3.5, 80, 50]*1e-4;
PolA = zeros(1,21); PolA(4) = [4.9]*1e-4;
Magnet{end}.Random{2} = struct( ...
    'PolynomB', PolB, ...     % All values are RMS, and specify variation of and relative to the main component
    'PolynomA', PolA, ...     % All values are RMS, and specify variation of and relative to the main component
    'Scaling', 1);          % Multiplicative factor, useful for current errors

PolB = zeros(1,21); PolB([9 15 21]) = 0.5e-4;
PolA = zeros(1,21);
Magnet{end}.Systematic{2} = struct( ...
    'PolynomB', PolB, ...     % All values are RMS
    'PolynomA', PolA, ...     % All values are RMS
    'Scaling', 1);          % Multiplicative factor, useful for current errors



% GENERIC OCTUPOLE MODEL
Magnet{end+1}.ID = 'Multipole';    % May be either atclass or FamName. Latter takes precedence if available.

Magnet{end}.Systematic{1}   =   DDR_ChallengingMachiningTolerances.Systematic;
Magnet{end}.Random{1}       =   DDR_ChallengingMachiningTolerances.Random;

Magnet{end}.Systematic{2} = struct( ...
    'PolynomB', [], ...   % All values are RMS
    'PolynomA', [], ...     % All values are RMS
    'Scaling', 1);          % Multiplicative factor, useful for current errors
Magnet{end}.Random{2} =     struct( ...
    'PolynomB', [], ...   % All values are RMS
    'PolynomA', [], ...     % All values are RMS
    'Scaling', 1);          % Multiplicative factor, useful for current errors



% GENERIC CORRECTOR MODEL
Magnet{end+1}.ID = 'Corrector';    % May be either atclass or FamName. Latter takes precedence if available.

Magnet{end}.Random{1} = DDR_ChallengingMachiningTolerances;



% SPECIFIC FAMILY ERRORS
Magnet{end+1}.ID = {'dipm','someBendyMag'};    % 'dipm' and 'someBendyMag', if they exist, will get the errors defined here rather than the generic ones for atclass = 'Bend'. 

Magnet{1}.Systematic{1}     =   DDR_ChallengingMachiningTolerances.Systematic;
Magnet{1}.Random{1}         =   DDR_ChallengingMachiningTolerances.Random;

% Random errors, incl. multipoles
PolB(1,[2 3 4 6 10]) = [2.5, 2.8, 1.9, 1.3, 0.3]*1e-4;
PolA(1,[3,4]) = [2.9, 1.4]*1e-4;
Magnet{end}.Random{2} = struct( ...
    'PolynomB', PolB, ...     % All values are RMS, and specify variation of and relative to the main component
    'PolynomA', PolA, ...     % All values are RMS, and specify variation of and relative to the main component
    'Scaling', .997);         % Multiplicative factor, useful for current errors

% Systematic errors, incl. multipoles
PolB(1,[6 10 14]) = [0.5, 0.5, 0.1]*1e-4;
PolA = [0,0,0,0];
Magnet{end}.Systematic{2} = struct( ...
    'PolynomB', PolB, ...     % All values are RMS
    'PolynomA', PolA, ...     % All values are RMS
    'Scaling', 1);            % Multiplicative factor, useful for current errors





% GENERIC BPM MODEL
Magnet{end+1}.ID = 'BPM';    % May be either atclass or FamName. Latter takes precedence if available.

Magnet{end}.Random{1}         =   DDR_ChallingingBPMCalibrationAccuracy.Random;


%% COLLECT THE OUTPUT
EM.Girder = Girder;
EM.Magnet = Magnet;

end

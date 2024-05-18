function EM = errormodel_DDRchallenging(varargin)
% Generates error model structure for use with "applyErrorModel"
%% Inputs 
% Optional arguments these are either scaling factors that multiply the 
%         errors defined in the challenging model defined below or (in the
%         case of strran) a relative excitation error
%
% gdran   = girder random alignment scaling factor, default = 1.0.
% mgalran = individual magnet random alignment scaling factor, defalt = 1.0.
% mulsys  = systematic magnet multipole error scaling factor, default = 1.0
% mulran  = random multipole error scaling factor, default = 1.0
% bpmran  = random bpm alignment error scaling factor, default = 1.0
% strran  = random magnet strength relative error, default = 0.02%
%% Usage examples
% ErrorModel=errormodel_DDRchallenging('gdran',0.0,'mgalran',1.0,'mulsys',0.0,'mulran',0.0,'strran',2E-4,'bpmran',0.0);
%

%% History
% 2024/05/16 : added input arguments as scalingf factors for each type of
%              error
%% Input argument parsing

gdran      = getoption(varargin,'gdran',1);
mgalran    = getoption(varargin,'mgalran',1);
mulsys     = getoption(varargin,'mulsys',1);
mulran     = getoption(varargin,'mulran',1);
bpmran     = getoption(varargin,'bpmran',1);
strran     = getoption(varargin,'strran',2E-4);



%% Girder error definitions
% Error model used is first based on girder name, which will be a


% Generic girder errors
Girder{1}.ID = {'Baseline','M1','M2','U1','U2','U3','U4','U5'};

Girder{1}.Systematic{1} = struct( ...
    'Pitch', 0, ...
    'Yaw', 0, ...
    'Roll', 0, ...
    'Sway', 0, ...
    'Heave', 0, ...
    'Surge', 0);
Girder{1}.Random{1} =     struct( ...
    'Pitch', 0e-6*gdran, ...
    'Yaw', 0e-6*gdran, ...
    'Roll', 0.1e-3*gdran, ...
    'Sway', 50e-6*gdran, ...
    'Heave', 50e-6*gdran, ...
    'Surge', 0*gdran);



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
    'Pitch', 0*mgalran, ...
    'Yaw', 0*mgalran, ...
    'Roll', 0.1e-3*mgalran, ...
    'Sway', 10e-6*mgalran, ...
    'Heave', 10e-6*mgalran, ...
    'Surge', 0*mgalran);
DDR_ChallingingBPMCalibrationAccuracy.Random = struct(...
    'Roll', 0.1e-3*bpmran, ...
    'Sway', 10e-6*bpmran, ...
    'Heave', 10e-6*bpmran);

% BASELINE FIELD ERROR, MAIN COMPONENT
StandardFieldError = strran;   % DDR assumes a field error of 0.02% RMS across the board, truncated at 2 sigma, for the main components


% GENERIC QUADRUPOLE MODEL
Magnet{1}.ID = 'Quadrupole';    % May be either atclass or FamName. Latter takes precedence if available.

Magnet{1}.Systematic{1}     =   DDR_ChallengingMachiningTolerances.Systematic;
Magnet{1}.Random{1}         =   DDR_ChallengingMachiningTolerances.Random;

% Random errors, incl. multipoles
PolB(1,[2 3 4 6 10]) = [2.5, 2.8, 1.9, 1.3, 0.3]*1e-4*mulran;
PolA(1,[3,4]) = [2.9, 1.4]*1e-4*mulran;
Magnet{end}.Random{2} = struct( ...
    'PolynomB', PolB, ...     % All values are RMS, and specify variation of and relative to the main component
    'PolynomA', PolA, ...     % All values are RMS, and specify variation of and relative to the main component
    'Scaling', 1 + StandardFieldError);            % Multiplicative factor, useful for current errors

% Systematic errors, incl. multipoles
PolB(1,[6 10 14]) = [0.5, 0.5, 0.1]*1e-4*mulsys;
PolA = [0,0,0,0]*mulsys;
Magnet{end}.Systematic{2} = struct( ...
    'PolynomB', PolB, ...     % All values are RMS
    'PolynomA', PolA, ...     % All values are RMS
    'Scaling', 1);            % Multiplicative factor, useful for current errors



% GENERIC SEXTUPOLE MODEL
Magnet{end+1}.ID = 'Sextupole';    % May be either atclass or FamName. Latter takes precedence if available.

Magnet{end}.Systematic{1}   =   DDR_ChallengingMachiningTolerances.Systematic;
Magnet{end}.Random{1}       =   DDR_ChallengingMachiningTolerances.Random;

PolB = zeros(1,21); PolB([3 4 5 9 15]) = [5.0, 5.2, 3.5, 80, 50]*1e-4*mulran;
PolA = zeros(1,21); PolA(4) = [4.9]*1e-4*mulran;
Magnet{end}.Random{2} = struct( ...
    'PolynomB', PolB, ...     % All values are RMS, and specify variation of and relative to the main component
    'PolynomA', PolA, ...     % All values are RMS, and specify variation of and relative to the main component
    'Scaling', 1 + StandardFieldError);          % Multiplicative factor, useful for current errors

PolB = zeros(1,21); PolB([9 15 21]) = 0.5e-4*mulsys;
PolA = zeros(1,21)*mulsys;
Magnet{end}.Systematic{2} = struct( ...
    'PolynomB', PolB, ...     % All values are RMS
    'PolynomA', PolA, ...     % All values are RMS
    'Scaling', 1);          % Multiplicative factor, useful for current errors



% GENERIC OCTUPOLE MODEL
Magnet{end+1}.ID = 'Multipole';    % May be either atclass or FamName. Latter takes precedence if available.

Magnet{end}.Systematic{1}   =   DDR_ChallengingMachiningTolerances.Systematic;
Magnet{end}.Random{1}       =   DDR_ChallengingMachiningTolerances.Random;

Magnet{end}.Random{2} =     struct( ...
    'PolynomB', []*mulran, ...   % All values are RMS
    'PolynomA', []*mulran, ...     % All values are RMS
    'Scaling', 1  + StandardFieldError);          % Multiplicative factor, useful for current errors

Magnet{end}.Systematic{2} = struct( ...
    'PolynomB', []*mulsys, ...   % All values are RMS
    'PolynomA', []*mulsys, ...     % All values are RMS
    'Scaling', 1);          % Multiplicative factor, useful for current errors


% GENERIC DIPOLE MODEL
Magnet{end+1}.ID = 'Bend';    % May be either atclass or FamName. Latter takes precedence if available.

% Generic machining misalignments are used
Magnet{end}.Systematic{1}   =   DDR_ChallengingMachiningTolerances.Systematic;
Magnet{end}.Random{1}       =   DDR_ChallengingMachiningTolerances.Random;

% PolB = zeros(1,21); PolB([3 4 5 9 15]) = [5.0, 5.2, 3.5, 80, 50]*1e-4;
% PolA = zeros(1,21); PolA(4) = [4.9]*1e-4;
Magnet{end}.Random{2} = struct( ...
    'PolynomB', []*mulran, ...     % All values are RMS, and specify variation of and relative to the main component
    'PolynomA', []*mulran, ...     % All values are RMS, and specify variation of and relative to the main component
    'Scaling', 1 + StandardFieldError);          % Multiplicative factor, useful for current errors

% PolB = zeros(1,21); PolB([9 15 21]) = 0.5e-4;
% PolA = zeros(1,21);
Magnet{end}.Systematic{2} = struct( ...
    'PolynomB', []*mulsys, ...     % All values are RMS
    'PolynomA', []*mulsys, ...     % All values are RMS
    'Scaling', 1);          % Multiplicative factor, useful for current errors



% GENERIC CORRECTOR MODEL
Magnet{end+1}.ID = 'Corrector';    % May be either atclass or FamName. Latter takes precedence if available.

Magnet{end}.Random{1} = DDR_ChallengingMachiningTolerances;



% GENERIC BPM MODEL
Magnet{end+1}.ID = 'BPM';    % May be either atclass or FamName. Latter takes precedence if available.

Magnet{end}.Random{1}         =   DDR_ChallingingBPMCalibrationAccuracy.Random;


%% COLLECT THE OUTPUT
EM.Girder = Girder;
EM.Magnet = Magnet;

end

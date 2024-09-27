function EM = errormodel_standard(varargin)

% Generates error model structure for use with "applyErrorModel"
%% Inputs 
% Optional arguments these are either scaling factors that multiply the 
%         errors defined in the chalenging model defined below or (in the
%         case of strran) a relative execitation error
%
% gdran   = girder random alignment scaling factor, default = 1.0.
% mgalran = individual magnet random alignment scaling factor, defalt = 1.0.
% mulsys  = systematic magnet multipole error scaling factor, default = 1.0
% mulran  = random multipole error scaling factor, default = 1.0
% bpmran  = random bpm alignment error scaling factor, default = 1.0
% strran  = random magnet strength relative error, default = 0.0

%% Usage examples
% ErrorModel=errormodel_standard('gdran',0.0,'mgalran',1.0,'mulsys',0.0,'mulran',0.0,'strran',0.0,'bpmran',0.0);
%

%% The section below seems to be an older code version that was later superseded
%  my suggestion is to turn this sectin into a documentation of what the 
%  "challenging" model is
% %% DEFINE ERRORS
%
% % ----------------------------------------
% % single magnet error table (RMS) --> MAGe
% % ----------------------------------------
% %      grad(frac)   dx(um)    dy(um)
% gradZero  = 1;  % 1 turn off/on the gradient errors
% shiftZero = 1;  % 1 turn off/on the displacement errors
% rollZero  = 1;  % 1 turn off/on the roll errors
% eQ = [ 1e-3*gradZero       20e-6*shiftZero   20e-6*shiftZero 20e-6*rollZero ]; % quadrupole
% eR = [ 1e-3*gradZero       20e-6*shiftZero   20e-6*shiftZero 20e-6*rollZero ]; % reverse-bends
% eS = [ 1e-3*gradZero       20e-6*shiftZero   20e-6*shiftZero 20e-6*rollZero ]; % sextupole
% eO = [ 1e-3*gradZero       20e-6*shiftZero   20e-6*shiftZero 20e-6*rollZero ]; % octupole
% eD = [ 1.2e-3*gradZero       20e-6*shiftZero   20e-6*shiftZero 20e-6*rollZero ]; % dipole
% MAGe.eQ = eQ; MAGe.eR = eR; MAGe.eS = eS; MAGe.eO = eO; MAGe.eD = eD;
% 
% % errors.Quadrupole =
% % errors.ReverseBend =
% % errors.Sextupole =
% % errors.Octupole =
% % errors.Bend = [ 0.07 0.16  1.2e-3];
% % errors.
% 
% % -------------------------------
% % girder random error table (RMS)
% % -------------------------------
% % <MSj> Given the girder shape, i.e. it's longer than it's wide, the
% % expectation is that roll will be harder to correctly determine. I
% % therefore swapped the yaw/pitch and roll values.
% 
% Girder.Pitch = [0 0 0 0 0 0 0];
% Girder.Yaw = [0 0 0 0 0 0 0];
% Girder.Roll = ones(1,7)*100e-6;
% Girder.Sway = ones(1,7)*50e-6;
% Girder.Heave = ones(1,7)*50e-6;
% Girder.Surge = zeros(1,7);
% 
% Magnet.Pitch = 0;
% Magnet.Yaw = 0;
% Magnet.Roll = 100e-6;
% Magnet.Sway = 25e-6;
% Magnet.Heave = 25e-6;
% Magnet.Surge = 0;
% 
% 
% Girder.Name
% Girder.Systematic
% Girder.Random
% 
% Magnet.Systematic
% Magnet.Random
% 
% What -- Filter -- Error class  -- Systematic/Random
% Girder  FamName   Field
% Magnet  Type      Alignment
% 
% 
% %        sway(um) heave(um) yaw(urad) pitch(urad) roll(urad)
% grdZero = 0.5;
% eGr{1}  = [100      100       10        10          25] * 1e-6  *grdZero;
% eGr{2}  = [100      100       10        10          25] * 1e-6  *grdZero;
% eGr{3}  = [100      100       10        10          25] * 1e-6  *grdZero;
% eGr{4}  = [100      100       10        10          25] * 1e-6  *grdZero;
% eGr{5}  = [100      100       10        10          25] * 1e-6  *grdZero;
% eGr{6}  = [100      100       10        10          25] * 1e-6  *grdZero;
% eGr{7}  = [100      100       10        10          25] * 1e-6  *grdZero;
% 
%
%% History
% 2024/03/29 : added input arguments as scalingf factors for each type of
%              error
%% Input argument parsing

gdran      = getoption(varargin,'gdran',1);
mgalran    = getoption(varargin,'mgalran',1);
mulsys     = getoption(varargin,'mulsys',1);
mulran     = getoption(varargin,'mulran',1);
bpmran     = getoption(varargin,'bpmran',1);
strran     = getoption(varargin,'strran',1);

%% Girder error definitions
% Error model used is first based on girder name, which will be a

% Generic girder errors
Girder{1}.ID = 'Baseline';

Girder{1}.Systematic{1} = struct('Pitch', 0, ...
    'Yaw', 0, ...
    'Roll', 0, ...
    'Sway', 0, ...
    'Heave', 0, ...
    'Surge', 0);
Girder{1}.Random{1} =     struct('Pitch', 0e-6*gdran, ...
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
DDR_ChallengingMachiningTolerances.Systematic = struct(...
    'Pitch', 0, ...
    'Yaw', 0, ...
    'Roll', 0, ...
    'Sway', 0, ...
    'Heave', 0, ...
    'Surge', 0);
DDR_ChallengingMachiningTolerances.Random = struct(...
    'Pitch', 0*mgalran , ...
    'Yaw', 0*mgalran , ...
    'Roll', 0.1e-3*mgalran, ...
    'Sway', 10e-6*mgalran, ...
    'Heave', 10e-6*mgalran, ...
    'Surge', 0*mgalran);
DDR_ChallingingBPMCalibrationAccuracy.Random = struct(...
    'Roll', 0.1e-3*bpmran, ...
    'Sway', 10e-6*bpmran, ...
    'Heave', 10e-6*bpmran);


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
    'Scaling', 1+strran);      % Multiplicative factor, useful for current errors

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
    'Scaling', 1+strran);          % Multiplicative factor, useful for current errors

PolB = zeros(1,21); PolB([9 15 21]) = 0.5e-4*mulsys;
PolA = zeros(1,21);
Magnet{end}.Systematic{2} = struct( ...
    'PolynomB', PolB, ...     % All values are RMS
    'PolynomA', PolA, ...     % All values are RMS
    'Scaling', 1);            % Multiplicative factor, useful for current errors


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
    'Scaling', 1+strran);          % Multiplicative factor, useful for current errors


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
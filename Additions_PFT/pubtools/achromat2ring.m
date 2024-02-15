function RING = achromat2ring(ACHR)
%ACHROMAT2RING expands an AT2 lattice to a full ring (360 degrees)
%
%[RING] = achromat2ring(ACHR)
%
%  INPUTS:
%  1. ACHROMAT  Initial AT structure
%
%  OUPUTS:
%  1. RING      Output ring
%
%  NOTES:
%  1. The function attempts to rely on the input achromat as far as
%  possible when determining global parameters. Fallback values will be the
%  MAX IV 3 GeV ring values. In particular, ensure the input achromat has   
%  total voltage and energy set.
%  2. Cavities and ringparam elements will not be duplicated
%  3. If no cavity exists, a thin cavity will be inserted at the end of the
%  lattice. It will attempt to get information from the input achromat if
%  available, i.e. from the RingParam element.
%  4. If ring parameters (circumference, periodicity, harmonic number)
%  match a known accelerator an attempt will be made to distribute extra
%  cavities in the correct achromats.
%
%  See also atenable_6d, atdisable_6d, atSetRingProperties



%% Default paramaters
E0 = 3e9;                   % The energy is assumed to be 3 GeV.
totalVoltage = 6*300e3;     % The existing MAX IV RF system is assumed to be kept
rfFreqGuess = 100e6;        % The existing MAX IV RF system is assumed to be kept


%% Update defaults if the lattice actually has the information
% Check if there's a mention of the frequency somewhere in the lattice.
% Note that the frequency will be recalculated regardless to match the
% circumference and energy, but we need to know whether it's a 100 MHz, 500
% MHz, or other frequency system. Else the harmonic number cannot be
% determined.
achrProps = atGetRingProperties(ACHR,'all');

% Particle energy. Check ringparam first, if that does not exist find the
% first mention of the energy inside the lattice.
if isfield(achrProps,'Energy') && ~isnan(achrProps.Energy), E0 = achrProps.Energy;
else
    I = findcells(ACHR,'Energy');
    if ~isempty(I)
        E0 = ACHR{I(1)}.Energy;
    end
end

% Frequency. Check ringparam first, otherwise scan the ring and grab the
% lowest frequency found, which hopefully corresponds to the main cavities.
if isfield(achrProps,'rf_frequency') && ~isnan(achrProps.rf_frequency), rfFreqGuess = achrProps.rf_frequency;
elseif ~isempty(findcells(ACHR,'Frequency'))
    I = findcells(ACHR,'Frequency');
    freqs = getcellstruct(ACHR,'Frequency',I);
    rfFreqGuess = min(freqs);
end

% Total voltage
if isfield(achrProps,'rf_voltage') && ~isnan(achrProps.rf_voltage), totalVoltage = achrProps.rf_voltage;
else
    I = findcells(ACHR,'Voltage');
    if ~isempty(I)
        totalVoltage = sum(getcellstruct(ACHR,'Voltage',I));
    end
end

%% Calculate and check some basics
% Note that the below assumes the electron actually has mass, which may be
% a problem. It's reset further down via a call to atSetRingProperties.
cellLength = findspos(ACHR,numel(ACHR)+1);
cellNumber = round(2*pi/sum(getcellstruct(ACHR,'BendingAngle',findcells(ACHR,'BendingAngle'))));
circumference = cellLength*cellNumber;

RelGamma    = E0/(PhysConstant.electron_mass_energy_equivalent_in_MeV.value * 1e6);
RelBeta     = sqrt(1 - RelGamma^-2);

harmNumber = round(circumference * rfFreqGuess/3e8);
cavFreq = harmNumber * RelBeta * physconst('LightSpeed')/circumference;

IS_6D=check_6d(ACHR);

%% Expand to a full ring

% Determine whether the input lattice contains a cavity, which is needed for some
% of the parameter calculations. If not, add a thin one. NB! Possible issue
% due to no ThinCavityPass.
cavpts = findcells(ACHR,'Frequency');
if isempty(cavpts)
    cavpts = numel(ACHR)+1;
    ACHR{cavpts} = atrfcavity('CAVITY',0,totalVoltage,cavFreq, harmNumber, E0, 'DriftPass');
end

% Duplicate it
N = ones(cellNumber,1)*(1:numel(ACHR)); N = N'; N = N(:);
RING = ACHR(N);

% MACHINE-DEPENDENT SECTION
% Default achromats using cavities are set here. A judgement is done based
% on circumference, periodicity (cellNumber) and harmonic number. With a
% few exceptions this should be a sufficient finger print (one example
% of an exception would be SOLARIS and MAX IV R1). 
if round(circumference) == 528 && cellNumber == 20 && harmNumber == 176
    cI = [1, 16:20];
end

% Get rid of any duplicated cavities by replacing them with drifts; we
% don't have cavities in all achromats
I = findcells(RING,'Frequency');
if exist('cI','var') && numel(I) == 20
        I = I(setdiff(1:20,cI));
else
    I = setdiff(I,findcells(ACHR,'Frequency'));
end
removeCavities(I);

    function removeCavities(I)
        for n = numel(I):-1:1
            if isfield(RING{I(n)},'Frequency')
                if RING{I(n)}.Length > 0
                    RING{I(n)} = atdrift('CavDriftReplacement',RING{I(n)}.Length);
                else
                    RING(I(n)) = [];
                end
            end
        end
    end

% Get rid of any duplicated RingParam elements, one is enough. Also, the
% new periodicity is 1.
I = findcells(RING,'Class','RingParam');
if ~isempty(I)
    RING = RING(setdiff(1:numel(RING),I(2:end)));
    RING{1}.Periodicity = 1;
end

% Ensure the ring properties are correctly set
RING = atSetRingProperties(RING,'HarmNumber',harmNumber,'particle','relativistic');     % This means relBeta == 0, which is implied to be an underlying assumption in AT tracking.
RING = atSetRingProperties(RING,'cavpts',findcells(RING,'Frequency'));
RING = atSetRingProperties(RING,'rf_frequency','nominal');                              % This does the usual frequency calculation, including relBeta. Which will be 1 now, unlike above.

% Set the new ring to have the same state as the incoming achromat
[~,RING]=check_6d(RING,IS_6D,'force');

end

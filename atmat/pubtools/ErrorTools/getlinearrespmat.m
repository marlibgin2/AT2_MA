function ORM = getlinearrespmat(varargin)
%GETLINEARRESPMAT Fast calculation of the orbit response matrix using LOCO
% ORM = getlinearrespmat(RING)

if nargin < 1
    error('getlinearrespmat:No lattice in input!');
end

if nargin > 0
    RINGData.Lattice = varargin{1};
end

if nargin > 1
    BPMData.BPMIndex = varargin{2};
else
    BPMData.BPMIndex = findcells(RINGData.Lattice,'Class','Monitor');
end

if nargin > 2
    CMData.HCMIndex = varargin{3};
else
    CMData.HCMIndex = findcells(RINGData.Lattice,'FamName','corrh');
end

if nargin > 2
    CMData.VCMIndex = varargin{4};
else
    CMData.VCMIndex = findcells(RINGData.Lattice,'FamName','corrv');
end

%% Default values
default_kick = 10e-6;
default_dpp = 1e-3;
alpha = mcf(RINGData.Lattice);

% Calculate frequency shift for the dispersive term
f0 = atGetRingProperties(RINGData.Lattice,'rf_frequency');
df = f0*-alpha*default_dpp;

RINGData.Lattice = RINGData.Lattice; 
RINGData.CavityFrequency = atGetRingProperties(RINGData.Lattice,'rf_frequency'); 
RINGData.CavityHarmNumber = atGetRingProperties(RINGData.Lattice,'harmonic_number');
CMData.HCMKicks = default_kick; 
CMData.VCMKicks = default_kick;
CMData.HCMCoupling = zeros(size(CMData.HCMIndex));
CMData.VCMCoupling = zeros(size(CMData.VCMIndex));

RM = locoresponsematrix(RINGData,BPMData,CMData,'rf',df);

RM = mat2cell(RM, ...
    [ numel(BPMData.BPMIndex), numel(BPMData.BPMIndex)], ...
    [ numel(CMData.HCMIndex), numel(CMData.VCMIndex), 1]);

% Translate to a format usable by atcorrectorbit
ORM.OrbHCor{1} = RM{1,1}./default_kick;
ORM.OrbHCor{3} = RM{2,1}./default_kick;
ORM.OrbVCor{1} = RM{1,2}./default_kick;
ORM.OrbVCor{3} = RM{2,2}./default_kick;
ORM.OrbHDPP = RM{1,3}./(df / (f0*-alpha));
ORM.OrbVDPP = RM{2,3}./(df / (f0*-alpha));
end

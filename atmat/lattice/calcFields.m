function FG = calcFields(varargin)
% Calculates the field and gradients of all the dipoles, quadrupoles,
% sextupoles and octupoles in T, T/m, T/m² and T/m³ respectively
% Dipole:     θ = BL/Bρ
% Quadrupole: K=(1/Bρ)*dB/dx
% sextupole:  S=(1/2)*(1/Bρ)*d²B/dx²
% Octupole:   O=(1/6)*(1/Bρ)*d³B/dx³
% 
%% Inputs
% Mandatory arguments
%   RING: AT2 lattice array (may be an achromat or full ring)
%   All_fams: cell array of string with the names of families for which
%             fields and gradients are to be determined
%
% Optional arguments
%   desc    : string describing the lattice, defaule = date
%   split   : factor by which to split the lattice elements, default = 1
%   verbose :defines level of verbose output, default=0, i.e. no output 
%
%% Outputs
% FG structure with the followiong fields 
%   FG.desc  : lattice description
%   FG.fams  : Name of all families given as inputs
%   FG.Spos  : Longitudinal position of the centre of all elements [m]
%
%   FG.Field : Dipole field B in [T] for all magnetic elements 
%   FG.GradQ : Quadrupole gradient dB/dx [T/m] for magnetic all elements 
%   FG.GradS : Sextupole gradient (1/2)*d²B/dx² [T/m²] for all magnetic elements
%   FG.GradO : Octupole gradient  (1/6)*d³B/dx³ [T/m³] for all magnetic elements
%   FG.FAll  : Dipole field B in [T] for all elements 
%
%% Usage examples
% All_fams = {'Q1';'Q2';'R1'; 'D2';'D3';'D1';'Q3';...
% 'Q4'; 'S1';'S2';'S3';'S4';'S5';'O1';'O2'; 'O3'};
% FG=calcFields(ACHR,All_fams,'verbose',1);

%% History
% SJ  2024/06/26: first version
% PFT 2024/06/27: added documentation + derived particle rigidity
%                 from lattice cell array
% SJ  2024/06/27: modified to get output as structure
% PFT 2024/06/28: modfied calculation of longitudinal coordinate,
%                 modified definition of gradients to match the one used in the 
%                 MagnetStrengthLimits excel file and matlab structure.
%                 added calculation of multipoles for dipole magnets
%                 handled the case of non-existing PolynomB elements
%                 added lattice description field
%                 added possibility to split the lattice
% PFT 2024/07/24: added outut of magnetic fields over all elements
%

%% Input argument parsing
[RING,All_fams] = getargs(varargin,[],{});

desc          = getoption(varargin,'desc',['calcFields: ' sprintf('%s', datetime)]);
verboselevel  = getoption(varargin,'verbose',0);
split         = getoption(varargin,'split',1);


%% Calculates fields and gradients
B_rho=atGetRingProperties(RING,'BRho');

% Initialize structure fields
FG.desc = desc;
FG.fams = {};
FG.Spos = [];
FG.Field = [];
FG.GradQ = [];
FG.GradS = [];
FG.GradO = [];


if (split>1)
    if (verboselevel>0)
            fprintf('%s calcFields: splittting input lattice into %3d slices \n', datetime, split);
    end
    RING=splitlat(RING,split);
end
FG.FAll = zeros(numel(RING),1);

% Iterate through the lattice elements
for i = 1:length(RING)
    element = RING{i};
    Spos = (findspos(RING,i)+findspos(RING,i+1))/2;
    if (ismember(element.FamName, All_fams)&&isfield(element,'PolynomB'))
        FG.fams = [FG.fams; {element.FamName}];
        FG.Spos = [FG.Spos; Spos];
        PolynomB=element.PolynomB;
        if(numel(PolynomB)<4)
            PolynomB=padarray(PolynomB,[0,4-numel(PolynomB)],'post');
        end

        if strcmp(element.PassMethod, 'BndMPoleSymplectic4Pass')
          
            B_field = (element.BendingAngle*B_rho)/element.Length;
            QP_gradient = PolynomB(2)*B_rho;
            SP_gradient = PolynomB(3)*B_rho;
            OP_gradient = PolynomB(4)*B_rho;
            FG.Field = [FG.Field; B_field];
            FG.GradQ = [FG.GradQ; QP_gradient];
            FG.GradS = [FG.GradS; SP_gradient];
            FG.GradO = [FG.GradO; OP_gradient];

        elseif strcmp(element.PassMethod, 'StrMPoleSymplectic4Pass')
            
            QP_gradient = PolynomB(2)*B_rho;
            SP_gradient = PolynomB(3)*B_rho;
            OP_gradient = PolynomB(4)*B_rho;

            FG.Field = [FG.Field; 0];
            FG.GradQ = [FG.GradQ; QP_gradient];
            FG.GradS = [FG.GradS; SP_gradient];
            FG.GradO = [FG.GradO; OP_gradient];
        end
    end
    
    if (isfield(element,'BendingAngle'))
        FG.FAll(i)=(element.BendingAngle*B_rho)/element.Length;
    end
end



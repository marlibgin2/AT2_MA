function LaHF = lattice2LaHF(varargin)
%LATTICE2LAFH Generates a Lattice Hardware File based on an AT2 lattice
%
%  LaHF = lattice2LaHF(RING)
%  LaHF = lattice2LaHF(AccStruct)
%
%  INPUT
%  1. RING      -   {cell array of structs} AT2 lattice
%  2. AccStruct -   {struct} Struct containing at minimum the AT2 lattice
%                   of full ring or an achromat. Name and Description
%                   fields are used if available.
%
%  OUTPUT
%  1. LaHF      -   Data extracted from the lattice, for future use.
%
%  NOTES
%  1. The intended use for this function is to translate a lattice file
%  into the preferred format for various engineering groups, who will
%  require the input for their work. This will vary greatly between
%  facilities.
%  2. The default output format is the one preferred by the MAX IV
%  Laboratory Technical Division (vacuum group, Design Office, etc.), which
%  is an Excel file. That said this function establishes a number of MATLAB
%  tables which could be written to other data formats is desired, although
%  the user will have to adapt the code.



%% Defaults
outputFile = [];
inputFile = [];

%% Input handling
% Currently lattice cell arrays are accepted. To come, .mat files with a
% number of saved variables to ease batch processing (needs automated
% output naming).

for n = nargin:-1:1
    if iscell(varargin{n})
        if isLattice(varargin{n})
            ACHR = varargin{n};
        elseif all(cellfun(@isfile,varargin{n}))
            % If a list of filenames has been received, recursive calls
            % will be made to establish the output.
        end
    end
    if isstruct(varargin{n})
        % Find the field containing the lattice
        I = structfun(@isLattice,varargin{n});
        if sum(I) > 0
            % If a lattice was found, extract it and attempt to get
            % descriptions and names as well
            fList = fieldnames(varargin{n});
            I = find(I,1,'first');
            ACHR = varargin{n}.(fList{I});
            I = contains(fList,'name','IgnoreCase',true);
            if any(I), LaHF.Name = varargin{n}.(fList{I}); end
            I = contains(fList,'descr','IgnoreCase',true);
            if any(I), LaHF.Description = varargin{n}.(fList{I}); end
        end
        %
    end
    if ischar(varargin{n})
        switch lower(varargin{n})
            case 'outputfile'
                outputFile = varargin{n+1};
            case 'inputfile'
                inputFile = varargin{n+1};
        end
    end
end



% Set a timestamp for the file generation
LaHF.Timestamp = datestr(now);

% If no name was available, as a last attempt try to check what the input
% variable names are whether any of them match the naming convention
if ~isfield(LaHF,'Name')
    candidateName = inputFile;
    if checkCandidate(candidateName)
        LaHF.Name = inputFile;
    else
        for k = 1:numel(nargin)
            candidateName = inputname(k);
            if checkCandidate(candidateName)
                LaHF.Name = candidateName;
                break;
            end
        end
    end
end

% Verify whether the 



%% Get lattice coordinates
LaHF.LatticeCoordinates_NoSlices = lattice2LatCoord(ACHR,'Condense');
LaHF.LatticeCoordinates_Detailed = lattice2LatCoord(ACHR,'NoCondense');


%% Get magnet parameter summary
LaHF.MagParams = lattice2MagParams(ACHR);

%% Get dipole slice tables
Dipoles = getDipoles(ACHR);
for n = 1:numel(Dipoles)
    LaHF.(Dipoles(n).Name) = dipole2slicemodel(Dipoles(n));
end

%% Get ring parameters

% Expand to a full ring
RING = achromat2ring(ACHR);

% Calculate the ring parameters
RP = ringpara(RING);
[E,~] = atx(RING);  % Can't trust atx for emittance values if the input isn't stellar

LaHF.RingParameters.Qx = num2str(RP.tunes(1));
LaHF.RingParameters.Qy = num2str(RP.tunes(2));
LaHF.RingParameters.Jx = num2str(RP.dampingJ(1));
LaHF.RingParameters.Js = num2str(RP.dampingJ(3));
LaHF.RingParameters.Emittance = sprintf('%.5g pm rad',RP.emittx*1e12);
LaHF.RingParameters.EnergySpread = sprintf('%.3g %%',RP.sigma_E*1e2);
LaHF.RingParameters.MomentumCompaction = sprintf('%g',RP.alphac);
LaHF.RingParameters.Betax = E(1).beta(1);
LaHF.RingParameters.Betay = E(1).beta(2);
LaHF.Data.GlobalParameters = RP;
LaHF.Data.MachineFunctions = E;




%% Get orbit
ORIGINAL = max4_simple_AT2;
n = findcells(ORIGINAL,'FamName','AchrEnd');
ORIGINAL = ORIGINAL(1:n);
ORIGINAL{1}.Periodicity = 20;

LaHF.Orbit.Existing = lattice2orbit(ORIGINAL);
LaHF.Orbit.Upgraded = lattice2orbit(ACHR);

%% Write to file
writeLaHF(LaHF,outputFile);

end



function writeLaHF(LaHF,varargin)

% To get the output filename, in order of preference:
% a) Use explicitly given input, if not empty
% b) Generate based on LaHF 'Name' field
% c) Ask the user
filename = [];
if nargin > 1
    filename = varargin{1};
end

if isempty(filename)
    if isfield(LaHF,'Name')
        DirectoryName = pwd;
        FileName = LaHF.Name;
        filename = fullfile(DirectoryName,[FileName '.xlsx']);
    end
end

if isempty(filename)
    [FileName, DirectoryName] = uiputfile('*.xlsx', 'Select file to save');
    if FileName == 0
        fprintf('No file saved.\n');
        return;
    else
        filename = fullfile(DirectoryName,FileName);
    end
end

tabNames = setdiff(fieldnames(LaHF),{'Name','Description','Timestamp','GlobalParameters','MachineFunctions'});

for n = 1:numel(tabNames)
    % For LaHF data fields useful in batch analysis, do not generate a tab
    if any(strcmpi(tabNames{n},{'Data'})), continue; end

    % If available, always write the name and description to each sheet for
    % easy identification
    if isfield(LaHF,'Name'), writecell({LaHF.Name},filename,'Sheet',tabNames{n},'Range','C2'); end
    if isfield(LaHF,'Description'), writecell({LaHF.Description},filename,'Sheet',tabNames{n},'Range','C3'); end

    % Proceed to insert the relevant information in each tab
    switch tabNames{n}
        case 'Orbit'
            writecell({'Existing ring'},filename,'Sheet',tabNames{n},'Range','C5');
            writecell({'Upgraded ring'},filename,'Sheet',tabNames{n},'Range','F5');
            writetable(LaHF.Orbit.Existing,filename,'Sheet',tabNames{n},'Range','C8');
            writetable(LaHF.Orbit.Upgraded,filename,'Sheet',tabNames{n},'Range','F8');
        case {'LatticeCoordinates','LatticeCoordinates_Detailed','LatticeCoordinates_NoSlices'}
            writecell(LaHF.(tabNames{n}).Header,filename,'Sheet',tabNames{n},'Range','C8');
            writetable(LaHF.(tabNames{n}).Data,filename,'Sheet',tabNames{n},'Range','C9','WriteVariableNames',false);
        case 'RingParameters'
            fNames = fieldnames(LaHF.RingParameters);
            for k = 1:numel(fNames)
                writecell(fNames(k),filename,'Sheet',tabNames{n},'Range',sprintf('C%d',6+k));
                writecell({LaHF.RingParameters.(fNames{k})},filename,'Sheet',tabNames{n},'Range',sprintf('D%d',6+k));
            end
        case 'Data'
            % Do nothing, this is for processing in MATLAB only. This line
            % shouldn't even be reachable, ever.
        otherwise
            writetable(LaHF.(tabNames{n}),filename,'Sheet',tabNames{n},'Range','C8');
    end
end

end

function Orbit = lattice2orbit(ACHR)

% Achromat total deflection
achrDeflection = sum(getcellstruct(ACHR,'BendingAngle',findcells(ACHR,'BendingAngle')));

% Ignore thin elements
I = cellfun(@(x) x.Length ~= 0,ACHR);
ACHR = ACHR(I);

% Calculate the X-Y position data
POSDATA = atgeometry(ACHR,1:numel(ACHR)+1,'HAngle',achrDeflection/2);
x = cat(1,POSDATA.x);
y = cat(1,POSDATA.y);  % Note the rescaling

Orbit = table(x,y,'VariableNames',{'x [m]','y [m]'});
end

function sliceModelTable = dipole2slicemodel(dipole)


N = numel(dipole.BEND);
% SliceName, Length, Angle, K, rho
% sliceName = getcellstruct(dipole.BEND,'FamName',1:numel(dipole.BEND));
length = getcellstruct(dipole.BEND,'Length',1:numel(dipole.BEND));
angle = getcellstruct(dipole.BEND,'BendingAngle',1:numel(dipole.BEND))*180/pi;
K = getcellstruct(dipole.BEND,'PolynomB',1:numel(dipole.BEND),2);
S = getcellstruct(dipole.BEND,'PolynomB',1:numel(dipole.BEND),3);
O = getcellstruct(dipole.BEND,'PolynomB',1:numel(dipole.BEND),4);
rho = length./getcellstruct(dipole.BEND,'BendingAngle',1:numel(dipole.BEND));

sliceName = cellfun(@getSliceName,dipole.BEND,'UniformOutput',false);

% % Append slice number
% for n = 1:numel(sliceName)
%     sliceName{n} = append(sliceName{n},'_',num2str(n));
% end

    

% Add the summary row
sliceName{N+1} = 'SUM';
length(N+1) = sum(length(1:N));
angle(N+1) = sum(angle(1:N));
K(N+1) = sum(K(1:N).*length(1:N))/length(N+1);   % Average field
S(N+1) = sum(S(1:N).*length(1:N))/length(N+1);
O(N+1) = sum(O(1:N).*length(1:N))/length(N+1);
rho(N+1) = nan;

sliceModelTable = table(sliceName, length, angle, K, S, O, rho, ...
    'VariableNames',{'Slice','Length [m]','Angle [deg]','K [m^-2]','S [m^-3]','O [m^-4]','rho [m]'});
end

function Dipoles = getDipoles(ACHR)

% Ignore zero length elements. Needed in order to determine continuous
% blocks (want to avoid assumptions regarding nbr of slices and similar).
I = cellfun(@(x) x.Length ~= 0,ACHR);
ACHR = ACHR(I);

% Get all magnet elements with positive bending field, and determine start
% and end indices for each block.
I = cellfun(@(x) isfield(x,'BendingAngle') && x.BendingAngle > 0,ACHR);
startIndex = find([0; diff(I)] == 1);
endIndex = find([diff(I); 0] == -1);

% Throw away redundant dipoles so there's only one of each type left
DipoleNames = cellfun(@getFamName,ACHR(startIndex),'UniformOutput',false);
[DipoleNames,J, ~] = unique(DipoleNames);
startIndex = startIndex(J); endIndex = endIndex(J);

% Create the output struct
Dipoles = struct('Name',[],'BEND',[]);
for n = 1:numel(DipoleNames)
    Dipoles(n).Name = DipoleNames{n};
    Dipoles(n).BEND = ACHR(startIndex(n):endIndex(n));
end



end

function MagParams = lattice2MagParams(ACHR)

% Get all elements that does not have positive bends (which should be the main dipoles)
M = findcells(ACHR,'PolynomB');
D = find(cellfun(@(x) isfield(x,'BendingAngle') && x.BendingAngle > 0, ACHR));
J = setdiff(M,D);

% Extract polynom B (padding may be needed for lengths varying between
% elements)
PolB = getcellstruct(ACHR,'PolynomB',J);
MaxOrder = max(cellfun(@numel,PolB));
PolB = cellfun(@cellpad ,PolB,'UniformOutput',false);

    function y = cellpad(x)
        y = x;
        xlen = numel(x);
        y((xlen+1):MaxOrder) = 0;
    end

% Gather up the parameters
Params = cat(2, ...
    getcellstruct(ACHR,'Length',J), ...
    cell2mat(PolB), ...
    getcellstruct(ACHR,'BendingAngle',J));
Params(isnan(Params(:))) = 0;
[Params, K] = unique(Params,'rows');
Names = cellfun(@getFamName,ACHR(J(K)),'UniformOutput',false);

% Determine the magnet type
MagType = cellfun(@getType,ACHR(J(K)),'UniformOutput',false);

    function magType = getType(elem)
        if isfield(elem,'BendingAngle') && elem.BendingAngle < 0
            magType = 'ReverseBend';
            return;
        else
            order = find(elem.PolynomB,1,'first');
            switch order
                case 1
                    magType = 'Corrector';
                case 2
                    magType = 'Quadrupole';
                case 3
                    magType = 'Sextupole';
                case 4
                    magType = 'Octupole';
            end
        end
    end


%     function output = getBendingAngle(elem)
%         if isfield(elem,'BendingAngle') output = elem.BendingAngle; else output = 0; end
%     end
%
% BendingAngle = cellfun(@getBendingAngle,ACHR);
% BendingAngle = BendingAngle(J);

[Names, I] = sort(Names);
Params = Params(I,:);
MagType = MagType(I);
% BendingAngle = BendingAngle(I,:);

MagParams = table(Names,Params(:,1), MagType, Params(:,3),Params(:,4),Params(:,5),Params(:,6)*180/pi);
MagParams.Properties.VariableNames = {'Magnet','Length [m]','Type','k [m^-2]','s [m^-3]','o [m^-4]','Angle [deg]'};

end

function LatticeCoordinates = lattice2LatCoord(ACHR,varargin)

% Input handling
CondenseFlag = true;
if nargin > 1
    for n = 1:numel(varargin)
        switch lower(varargin{n})
            case 'condense'
                CondenseFlag = true;
            case 'nocondense'
                CondenseFlag = false;
        end
    end
end

% Get the particle energy
I = findcells(ACHR,'Energy');
E0 = unique(getcellstruct(ACHR,'Energy',I));
if numel(E0) > 1
    warning('lattice2LaHF:lattice2LatCoord:Several particle energies found in the lattice! Using the first one only.');
    E0 = E0(1);
end
Brho = -(E0 + PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e6)/PhysConstant.speed_of_light_in_vacuum.value;

% Ignore zero length elements
I = cellfun(@(x) x.Length ~= 0,ACHR);
ACHR = ACHR(I);

% Establish the output struct array
LatticeCoordinates.Data = struct('Name',[],'Start',[],'Centre',[],'End',[],'Length',[],'BendingAngle',[],'Field',[]);
LatticeCoordinates.Header = {'Name','Start [m]','Centre [m]','End [m]','Length [m]','Bending angle [deg]','By @ x=0 [T]'};

s = findspos(ACHR,1:numel(ACHR)+1);
theta = cellfun(@getTheta,ACHR);

    function y = getTheta(x)
        if isfield(x,'BendingAngle')
            y = x.BendingAngle * 180/pi;
        else
            y = 0;
        end
    end

% Name, start, end, centre, length
% Element Name, Start Position, Centre Position, End Position , Length, Bending Angle.

for n = 1:numel(ACHR)
    %     if ACHR{n}.Length == 0, ACHR{n} = []; continue; end

    switch ACHR{n}.PassMethod
        case 'DriftPass'
            LatticeCoordinates.Data(n).Name = 'Straight';
            LatticeCoordinates.Data(n).Start = s(n);
            LatticeCoordinates.Data(n).End = s(n+1);
            LatticeCoordinates.Data(n).Centre = (LatticeCoordinates.Data(n).Start + LatticeCoordinates.Data(n).End)/2;
            LatticeCoordinates.Data(n).Length = s(n+1)-s(n);
            LatticeCoordinates.Data(n).BendingAngle = 0;
            LatticeCoordinates.Data(n).Field = 0;
        case {'BendLinearPass','BndMPoleSymplectic4E2Pass','BndMPoleSymplectic4E2RadPass','BndMPoleSymplectic4Pass','BndMPoleSymplectic4QuantPass','BndMPoleSymplectic4RadPass','BndStrMPoleSymplectic4Pass'}
            if CondenseFlag, LatticeCoordinates.Data(n).Name = getFamName(ACHR{n});
            else, LatticeCoordinates.Data(n).Name = getSliceName(ACHR{n}); end
            LatticeCoordinates.Data(n).Start = s(n);
            LatticeCoordinates.Data(n).End = s(n+1);
            LatticeCoordinates.Data(n).Centre = (LatticeCoordinates.Data(n).Start + LatticeCoordinates.Data(n).End)/2;
            LatticeCoordinates.Data(n).Length = s(n+1)-s(n);
            LatticeCoordinates.Data(n).BendingAngle = theta(n);
            LatticeCoordinates.Data(n).Field = Brho * (theta(n)*pi/180) / (s(n+1)-s(n));
        case {'QuadLinearFPass','QuadLinearPass','QuantDiffPass','StrMPoleSymplectic4Pass','StrMPoleSymplectic4QuantPass','StrMPoleSymplectic4RadPass','ThinMPolePass','VariableThinMPolePass'}
            if CondenseFlag, LatticeCoordinates.Data(n).Name = getFamName(ACHR{n});
            else, LatticeCoordinates.Data(n).Name = getSliceName(ACHR{n}); end
            LatticeCoordinates.Data(n).Start = s(n);
            LatticeCoordinates.Data(n).End = s(n+1);
            LatticeCoordinates.Data(n).Centre = (LatticeCoordinates.Data(n).Start + LatticeCoordinates.Data(n).End)/2;
            LatticeCoordinates.Data(n).Length = s(n+1)-s(n);
            LatticeCoordinates.Data(n).BendingAngle = 0;
            LatticeCoordinates.Data(n).Field = 0;
            %     case {'IdentityPass','AperturePass'}
            %         continue;
        case     {'CorrectorPass'               }
            if CondenseFlag, LatticeCoordinates.Data(n).Name = getFamName(ACHR{n});
            else, LatticeCoordinates.Data(n).Name = getSliceName(ACHR{n}); end
            LatticeCoordinates.Data(n).Start = s(n);
            LatticeCoordinates.Data(n).End = s(n+1);
            LatticeCoordinates.Data(n).Centre = (LatticeCoordinates.Data(n).Start + LatticeCoordinates.Data(n).End)/2;
            LatticeCoordinates.Data(n).Length = s(n+1)-s(n);
            LatticeCoordinates.Data(n).BendingAngle = 0;
            LatticeCoordinates.Data(n).Field = 0;
        otherwise
            error('lattice2LaHF:unknown element type! Could not produce LaHF.');
    end

end

% If slices should be condensed, and only an equivalent bulk field
% reported...
if CondenseFlag
    for n = numel(LatticeCoordinates.Data):-1:2
        if strcmp(LatticeCoordinates.Data(n).Name,LatticeCoordinates.Data(n-1).Name)
            LatticeCoordinates.Data(n-1).End = LatticeCoordinates.Data(n).End;
            LatticeCoordinates.Data(n-1).Length = LatticeCoordinates.Data(n-1).End - LatticeCoordinates.Data(n-1).Start;
            LatticeCoordinates.Data(n-1).Centre = (LatticeCoordinates.Data(n-1).Start + LatticeCoordinates.Data(n-1).End)/2;
            if isfield(LatticeCoordinates.Data(n),'BendingAngle')
                LatticeCoordinates.Data(n-1).BendingAngle = LatticeCoordinates.Data(n-1).BendingAngle + LatticeCoordinates.Data(n).BendingAngle;
            end
            LatticeCoordinates.Data(n) = [];
            n = n-1;
        end
    end
    LatticeCoordinates.Data = rmfield(LatticeCoordinates.Data,'Field');
    LatticeCoordinates.Header = setdiff(LatticeCoordinates.Header, 'By @ x=0 [T]','stable');
end

LatticeCoordinates.Data = struct2table(LatticeCoordinates.Data);



%
%
%     {'BeamLoadingCavityPass'       }
%     {'BeamMomentsPass'             }
%     {'CavityPass'                  }
%     {'DeltaQPass'                  }
%     {'EAperturePass'               }
%     {'ExactHamiltonianPass'        }
%     {'GWigSymplecticPass'          }
%     {'GWigSymplecticRadPass'       }
%     {'IdTablePass'                 }
%
%     {'ImpedanceTablePass'          }
%     {'Matrix66Pass'                }
%     {'MatrixTijkPass'              }
%     {'RFCavityPass'                }
%     {'SolenoidLinearPass'          }
%     {'WakeFieldPass'               }
%     {'WiggLinearPass'              }
end


%% Common helper functions

function y = getFamName(x)
if isstruct(x) && isfield(x,'FamName')
        x = x.FamName;
end
y = regexp(x,'_','split');
y = y{1};
end

function y = getSliceName(x)
if isstruct(x)
    if isfield(x,'SliceName')
        x = x.SliceName; 
    elseif isfield(x,'FamName')
        x = x.FamName;
    end
end
y = regexp(x,'_','split');
y = y{1};
end

function flag = isLattice(var)
if iscell(var)
    flag = all(cellfun(@isstruct,var)) && numel(var) > 2;
else
    flag = false;
end
end

% Highly facility-specific file, which will require adaptation based on the
% naming convention used.
function okFlag = checkCandidate(candidateName)
if isempty(candidateName), okFlag = false; return; end
candidateName = regexpi(candidateName,'m4U.[0-9]{6,8}.[a-z]*[\S0-9]*','match');
okFlag = ~isempty(candidateName);
end
function LaHF = lattice2LaHF(ACHR)

%% Input handling
% Currently lattice cell arrays are accepted. To come, .mat files with a
% number of saved variables to ease batch processing (needs automated
% output naming).

%% Get lattice coordinates
LaHF.LatticeCoordinates = lattice2LatCoord(ACHR);

%% Get magnet parameter summary
LaHF.MagParams = lattice2MagParams(ACHR);

%% Get dipole slice tables
Dipoles = getDipoles(ACHR);
for n = 1:numel(Dipoles)
    LaHF.(Dipoles(n).Name) = dipole2slicemodel(Dipoles(n));
end

%% Get ring parameters
% On to-do list.


%% Get orbit
ORIGINAL = max4_simple_AT2;
n = findcells(ORIGINAL,'FamName','AchrEnd');
ORIGINAL = ORIGINAL(1:n);
ORIGINAL{1}.Periodicity = 20;

LaHF.Orbit.Existing = lattice2orbit(ORIGINAL);
LaHF.Orbit.Upgraded = lattice2orbit(ACHR);

%% Write to file
writeLaHF(LaHF);

end


function writeLaHF(LaHF,varargin)
    
filename = [];
if nargin > 1
    filename = varargin{1};
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

tabNames = fieldnames(LaHF);

for n = 1:numel(tabNames)
%             writecell({'Existing ring'},filename,'Sheet',tabNames{n},'Range','C5');
%             writecell({'Upgraded ring'},filename,'Sheet',tabNames{n},'Range','F5');

    switch tabNames{n}
        case 'Orbit'
            writecell({'Existing ring'},filename,'Sheet',tabNames{n},'Range','C5');
            writecell({'Upgraded ring'},filename,'Sheet',tabNames{n},'Range','F5');
            writetable(LaHF.Orbit.Existing,filename,'Sheet',tabNames{n},'Range','C8');
            writetable(LaHF.Orbit.Upgraded,filename,'Sheet',tabNames{n},'Range','F8');
        case 'LatticeCoordinates'
            writetable(LaHF.LatticeCoordinates,filename,'Sheet',tabNames{n},'Range','C8');
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
y = cat(1,POSDATA.y)*1000;  % Note the rescaling

Orbit = table(x,y,'VariableNames',{'x [m]','y [mm]'});
end

function sliceModelTable = dipole2slicemodel(dipole)


N = numel(dipole.BEND);
% SliceName, Length, Angle, K, rho
sliceName = getcellstruct(dipole.BEND,'FamName',1:numel(dipole.BEND));
length = getcellstruct(dipole.BEND,'Length',1:numel(dipole.BEND));
angle = getcellstruct(dipole.BEND,'BendingAngle',1:numel(dipole.BEND))*180/pi;
K = getcellstruct(dipole.BEND,'PolynomB',1:numel(dipole.BEND),2);
S = getcellstruct(dipole.BEND,'PolynomB',1:numel(dipole.BEND),3);
O = getcellstruct(dipole.BEND,'PolynomB',1:numel(dipole.BEND),4);
rho = length./getcellstruct(dipole.BEND,'BendingAngle',1:numel(dipole.BEND));

% Add the summary row
sliceName{N+1} = 'SUM';
length(N+1) = sum(length(1:N));
angle(N+1) = sum(angle(1:N));
K(N+1) = sum(K(1:N).*length(1:N))/length(N+1);   % Average field
S(N+1) = sum(S(1:N).*length(1:N))/length(N+1);
O(N+1) = sum(O(1:N).*length(1:N))/length(N+1);
rho(N+1) = nan;

sliceModelTable = table(sliceName, length, angle, K, S, O, rho, ...
    'VariableNames',{'Slice','Length [mm]','Angle [deg]','K [m^-2]','S [m^-3]','O [m^-4]','rho [m]'});
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
    getcellstruct(ACHR,'Length',J)*1000, ...
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
MagParams.Properties.VariableNames = {'Magnet','Length [mm]','Type','k [m^-2]','s [m^-3]','o [m^-4]','Angle [deg]'};

end

function LatticeCoordinates = lattice2LatCoord(ACHR)

% Ignore zero length elements
I = cellfun(@(x) x.Length ~= 0,ACHR);
ACHR = ACHR(I);

% Establish the output struct array
LatticeCoordinates = struct('Name',[],'Start',[],'Centre',[],'End',[],'Length',[]);

s = findspos(ACHR,1:numel(ACHR)+1);

% Name, start, end, centre, length


for n = 1:numel(ACHR)
    %     if ACHR{n}.Length == 0, ACHR{n} = []; continue; end

    switch ACHR{n}.PassMethod
        case 'DriftPass'
            LatticeCoordinates(n).Name = 'Straight';
            LatticeCoordinates(n).Start = s(n);
            LatticeCoordinates(n).End = s(n+1);
            LatticeCoordinates(n).Centre = (LatticeCoordinates(n).Start + LatticeCoordinates(n).End)/2;
            LatticeCoordinates(n).Length = s(n+1)-s(n);
        case {'BendLinearPass','BndMPoleSymplectic4E2Pass','BndMPoleSymplectic4E2RadPass','BndMPoleSymplectic4Pass','BndMPoleSymplectic4QuantPass','BndMPoleSymplectic4RadPass','BndStrMPoleSymplectic4Pass'}
            LatticeCoordinates(n).Name = ACHR{n}.FamName;
            LatticeCoordinates(n).Start = s(n);
            LatticeCoordinates(n).End = s(n+1);
            LatticeCoordinates(n).Centre = (LatticeCoordinates(n).Start + LatticeCoordinates(n).End)/2;
            LatticeCoordinates(n).Length = s(n+1)-s(n);
        case {'QuadLinearFPass','QuadLinearPass','QuantDiffPass','StrMPoleSymplectic4Pass','StrMPoleSymplectic4QuantPass','StrMPoleSymplectic4RadPass','ThinMPolePass','VariableThinMPolePass'}
            LatticeCoordinates(n).Name = ACHR{n}.FamName;
            LatticeCoordinates(n).Start = s(n);
            LatticeCoordinates(n).End = s(n+1);
            LatticeCoordinates(n).Centre = (LatticeCoordinates(n).Start + LatticeCoordinates(n).End)/2;
            LatticeCoordinates(n).Length = s(n+1)-s(n);
            %     case {'IdentityPass','AperturePass'}
            %         continue;
        case     {'CorrectorPass'               }
            LatticeCoordinates(n).Name = ACHR{n}.FamName;
            LatticeCoordinates(n).Start = s(n);
            LatticeCoordinates(n).End = s(n+1);
            LatticeCoordinates(n).Centre = (LatticeCoordinates(n).Start + LatticeCoordinates(n).End)/2;
            LatticeCoordinates(n).Length = s(n+1)-s(n);

        otherwise
            error('lattice2LaHF:unknown element type! Could not produce LaHF.');
    end

end

    for n = numel(LatticeCoordinates):-1:2
        if strcmp(LatticeCoordinates(n).Name,LatticeCoordinates(n-1).Name)
            LatticeCoordinates(n-1).End = LatticeCoordinates(n).End;
            LatticeCoordinates(n-1).Length = LatticeCoordinates(n-1).End - LatticeCoordinates(n-1).Start;
            LatticeCoordinates(n-1).Centre = (LatticeCoordinates(n-1).Start + LatticeCoordinates(n-1).End)/2;
            LatticeCoordinates(n) = [];
            n = n-1;
        end
    end

    LatticeCoordinates = struct2table(LatticeCoordinates);



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
y = regexp(x.FamName,'_','split');
y = y{1};
end




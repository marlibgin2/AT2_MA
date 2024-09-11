function Fams = findFams(LAT)
% generates a structure containing info on the various types oe 
% element families in an AT2.0 cell array
%% Inputs
% LAT : AT2.0 lattice cell array
%
%% Outputs
% Fams: structure with the fields
% Fams.Dipoles    : cell array of strings with names of dipole families
% Fams.Multipoles : cell array of strings with names of multipole families
% Fams.Correctors : cell array of strings with names of orbit corrector
%                   families
% Fams.BPMs       : cell array of strings with names of monitor families
%                   
% Fams.nelem.Dipoles : array with number of elements belonging to eachoneof
%                      the dipole families
% Fams.nelem.Multipoles : array with number of elements belonging to each one of
%                      the multipole families
% Fams.nelem.Correctos : array with number of elements belonging to each one of
%                      the corrector families
% Fams.nelem.Correctos : array with number of elements belonging to each one of
%                      the corrector families
% Fams.nelem.BPMs : array with number of elements belonging to each one of
%                      the monitor families
%
FamTypes={'Dipoles';'Multipoles';'Correctors';'BPMs'}; 
nfamtypes=numel(FamTypes);
k=zeros(nfamtypes,1);
for i=1:nfamtypes
    Fams.(FamTypes{i})={};
end

for i=1:numel(LAT)
    EL=LAT{i};
%     PassMethod = atgetfieldvalues(LAT, i, 'PassMethod');
    if isfield(EL,'PassMethod'), PassMethod = {EL.PassMethod}; else, PassMethod = {}; end
    switch PassMethod{1}
        case {'BndMPoleSymplectic4Pass','BndMPoleSymplectic4RadPass'}
            k(1)=k(1)+1;
            Fams.Dipoles{k(1),1}=EL.FamName;

        case {'StrMPoleSymplectic4Pass', 'StrMPoleSymplectic4RadPass','ThinMPolePass'}
            k(2)=k(2)+1;
            Fams.Multipoles{k(2),1}=EL.FamName;
            
        case 'CorrectorPass'
            k(3)=k(3)+1;
            Fams.Correctors{k(3),1}=EL.FamName;

        case 'IdentityPass'
            Class=EL.Class;
            if (strcmpi(Class,'Monitor'))
                k(4)=k(4)+1;
                Fams.BPMs{k(4),1}=EL.FamName;
            end
        otherwise
    end
end

for i=1:nfamtypes
    totList = Fams.(FamTypes{i});
    FamList = unique(totList);
    for j = 1:numel(FamList)
        NumElems = sum(strcmp(FamList{j},totList));
        [Fams.(FamTypes{i})] = FamList;
        Fams.nelems.(FamTypes{i})(j,1) = NumElems;
    end
end

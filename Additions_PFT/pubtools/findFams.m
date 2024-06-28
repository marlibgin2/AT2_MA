function Fams = findFams(LAT)
% generates a structure containing info on the various types oe 
% element families in an AT2.0 cell array
%% Inputs
% LAT : AT2.0 lattice cell array
%
%% Outputs
% Fams: struncture with the fields
% Fams.Dipoles    : cell array of strings with names of dipole families
% Fams.Multipoles : cell aray of stris with names of multipole families
% Fams.Correctors : cell array of strings with names of orbit cOrrector
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
    PassMethod = atgetfieldvalues(LAT, i, 'PassMethod');
    switch PassMethod{1}
        case 'BndMPoleSymplectic4Pass'
            k(1)=k(1)+1;
            Fams.Dipoles{k(1),1}=EL.FamName;

        case 'StrMPoleSymplectic4Pass'
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
    Fams.(FamTypes{i})=unique(Fams.(FamTypes{i}));
    k(i)=numel(Fams.(FamTypes{i}));
    Fams.nelems.(FamTypes{i})=zeros(k(i),1);
    for j=1:k(i)
        I_fam = find(atgetcells(LAT, 'FamName', Fams.(FamTypes{i}){j}));
        Fams.nelems.(FamTypes{i})(j)=numel(I_fam);
    end
end

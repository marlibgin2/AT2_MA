function Fams = findFams(LAT)
% generates a cell array of strings containign all magnet families in 
% a AT2 lattice
%
Fams.Dipoles={};
Fams.Multipoles={};
Fams.Correctors={};
Fams.BPMs={};
kd=0;
km=0;
kc=0;
kb=0;

for i=1:numel(LAT)
    EL=LAT{i};
    PassMethod = atgetfieldvalues(LAT, i, 'PassMethod');
    switch PassMethod{1}
        case 'BndMPoleSymplectic4Pass'
            kd=kd+1;
            Fams.Dipoles{kd}=EL.FamName;

        case 'StrMPoleSymplectic4Pass'
            km=km+1;
            Fams.Multipoles{km}=EL.FamName;
            
        case 'CorrectorPass'
            kc=kc+1;
            Fams.Correctors{kc}=EL.FamName;
        case 'IdentityPass'
            Class=EL.Class;
            if (strcmpi(Class,'Monitor'))
                kb=kb+1;
                Fams.BPMs{kb}=EL.FamName;
            end
        otherwise
    end
end
    
Fams.Dipoles=unique(Fams.Dipoles);
Fams.Multipoles=unique(Fams.Multipoles);
Fams.Correctors=unique(Fams.Correctors);
Fams.BPMs=unique(Fams.BPMs);
Fams.ndipoles=nd;
Fams.nmultipoles=nm;
Fams.ncorrectors=nc;
Fams.bpms=nb;
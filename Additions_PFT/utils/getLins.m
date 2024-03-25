function DVs = getLins(nLAT, LAT, LatticeOptData)
% Gets the decision variables for quadrupole gradients from a given lattice
% Uses info in structure LatticeOptData to define what these variables are
% Returns a one line matrix of nvars values
%
% nLAT is :
%          1 if DVs are to be taken from LAT with same structure as HACHRO in LatticeOptData
%          2 if DVs are to be taken from LAT with same structure as ACHRO in LatticeOptData
%          3 if DVs are to be taken from LAT with same structure as UC in LatticeOptData
%          4 if DVs are to be taken from LAT with same structure as IMC1 in LatticeOptData
%          5 if DVs are to be taken from LAT with same structure as RING in LatticeOptData
%          6 if DVS are to be taken from LAT with same structure as RINGGRD in LatticeOptData
%          7 if DVS are to be taken from LAT with same structure as ACHROGRD in LatticeOptData
%
% Note : the calling routine must check that the choice of nLAT and LAT are
% compatible with each other.

optMode = LatticeOptData.optMode;

%
% check for backward compatibility. 
%
if (isfield(LatticeOptData,'nvars_lin'))
    nvars = LatticeOptData.nvars_lin;
else
    nvars = 7;
end
    
if (isfield(LatticeOptData,'stdfamlist_lin'))
    stdfamlist  = LatticeOptData.stdfamlist_lin;
else
    switch optMode
        case 'SIMP'
            stdfamlist = 1:7;
        case {'COMP','CHRO'}
            stdfamlist = [1 2 3 6 7];
    end
end

if (isfield(LatticeOptData,'nstdfamlist_lin'))
    nstdfamlist = LatticeOptData.nstdfamlist_lin;
else
    switch  optMode
        case 'SIMP'
             nstdfamlist = [];
        case {'COMP','CHRO'}
            nstdfamlist = [4 5];
    end        
end

if(isfield(LatticeOptData,'famtype_lin'))
    famtype=LatticeOptData.famtype_lin;
else
    famtype=ones(1,nvars);
end


UC          = LatticeOptData.UC;
ACHRO       = LatticeOptData.ACHRO;
HACHRO      = LatticeOptData.HACHRO;
RING        = LatticeOptData.RING;

if (isfield(LatticeOptData,'IMC1'))
    IMC1  = LatticeOptData.IMC1;
end

if (isfield(LatticeOptData,'RINGGRD'))
    RINGGRD = LatticeOptData.RINGGRD;
end

if (isfield(LatticeOptData,'ACHROGRD'))
    ACHROGRD = LatticeOptData.ACHROGRD;
end

DVs = NaN(1, nvars);

switch nLAT
    case 1
        Ifams  = LatticeOptData.IfamsHlin;
        if (length(LAT)~=length(HACHRO))
            fprintf('Warning: Incompatible input to getLins for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nvars);
            return
        end

    case 2
        Ifams  = LatticeOptData.IfamsFlin;
        if (length(LAT)~=length(ACHRO))
            fprintf('Warning: Incompatible input to getLins for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nvars);
            return
        end

    case 3
        Ifams  = LatticeOptData.IfamsUClin;
        if (length(LAT)~=length(UC))
            fprintf('Warning: Incompatible input to getLins for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nvars);
            return
        end

    case 4
        Ifams  = LatticeOptData.IfamsIMC1lin;
        if (length(LAT)~=length(IMC1))
            fprintf('Warning: Incompatible input to getLins for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nvars);
            return
        end

     case 5
        Ifams  = LatticeOptData.IfamsRINGlin;
        if (length(LAT)~=length(RING))
            fprintf('Warning: Incompatible input to getLins for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nvars);
            return
        end

    case 6
        Ifams  = LatticeOptData.IfamsRINGGRDlin;
        if (length(LAT)~=length(RINGGRD))
            fprintf('Warning: Incompatible input to getLins for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nallfams);   
            return
        end

    case 7
        Ifams  = LatticeOptData.IfamsACHROGRDlin;
        if (length(LAT)~=length(ACHROGRD))
            fprintf('Warning: Incompatible input to getLins for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nvars);      
            return
        end

    otherwise
            fprintf('%s Warning: Error in getLins, unknow lattice type nLAt = %2d \n',datetime, nLAT);
         
end
    
for i=1:length(stdfamlist)
    if (not(isempty(Ifams{stdfamlist(i)})))
        DVs(stdfamlist(i))=LAT{Ifams{stdfamlist(i)}(1)}.PolynomB(1,famtype(stdfamlist(i))+1);
    else
        DVs(stdfamlist(i))= NaN;
    end
end

for i=1:length(nstdfamlist)
     Ks = NaN(1,size(Ifams{nstdfamlist(i)},2));
     for l=1:size(Ifams{nstdfamlist(i)},1)
            Ks(l)= LAT{Ifams{nstdfamlist(i)}(l)}.PolynomB(1,2);
     end
     
     [kmax,pos] = max(abs(Ks));
     DVs(nstdfamlist(i))=sign(Ks(pos))*kmax;
end

end


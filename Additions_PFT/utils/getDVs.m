function DVs = getDVs(nLAT, LAT, LatticeOptData)
% Gets the decision variables from a given lattice
% Uses info in strcuture LatticeOptData to define what these variables are
% Returns a one line matrix of nvars values
%
% nLAT is :
%          1 if DVs are to be taken from LAT with same structure as HACHRO in LatticeOptData
%          2 if DVs are to be taken from LAT with same structure as ACHRO in LatticeOptData
%          3 if DVs are to be taken from LAT with same structure as UC in LatticeOptData
%          4 if DVs are to be taken from LAT with same structure as IMC1 in LatticeOptData
%          5 if DVs are to be taken from LAT with same structure as RING in LatticeOptData
%          6 if DVs are to be taken from LAT with same structure as RINGGRD in LatticeOptData
%          7 if DVs are to be taken from LAT with same structure as ACHROGRD in LatticeOptData
%
% Note : the calling routine must check that the choice of nLAT and LAT are
% compatible with each other.

%% History
% PFT 2023,first version
% PFT 2024/08/27: included bend angles (with fixed profile) and distance as
%                 DV
% PFT 2024/08/28: included slice bend angles as DVs

%% Preamble
nvars       = LatticeOptData.nvars;
%
% check for backward compatibility
%
%
if (isfield(LatticeOptData,'stdfamlist'))
    stdfamlist  = LatticeOptData.stdfamlist;
else
    stdfamlist = 1:7;
end

if (isfield(LatticeOptData,'nstdfamlist'))
    nstdfamlist = LatticeOptData.nstdfamlist;
else
    nstdfamlist = [];
end

if(isfield(LatticeOptData,'famtype'))
    famtype=LatticeOptData.famtype;
else
    famtype=ones(1,nvars);
end

if (isfield(LatticeOptData,'bafamlist'))
    bafamlist  = LatticeOptData.bafamlist;
else
    bafamlist = [];
end

if (isfield(LatticeOptData,'Lfamlist'))
    Lfamlist  = LatticeOptData.Lfamlist;
else
    Lfamlist = [];
end

if (isfield(LatticeOptData,'slicefamlist'))
    slicefamlist = LatticeOptData.slicefamlist;
else
    slicefamlist = [];
end

UC          = LatticeOptData.UC;
ACHRO       = LatticeOptData.ACHRO;
HACHRO      = LatticeOptData.HACHRO;

if (isfield(LatticeOptData,'IMC1'))
    IMC1  = LatticeOptData.IMC1;
end

if (isfield(LatticeOptData,'RING'))
    RING = LatticeOptData.RING;
end

if (isfield(LatticeOptData,'RINGGRD'))
    RINGGRD = LatticeOptData.RINGGRD;
end

if (isfield(LatticeOptData,'ACHROGRD'))
    ACHROGRD = LatticeOptData.ACHROGRD;
end
DVs = NaN(1, nvars);

%% Finds lattice locations
switch nLAT
    case 1
        Ifams  = LatticeOptData.IfamsH;
        if (length(LAT)~=length(HACHRO))
            fprintf('Warning: Incompatible input to getDVs for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nvars);
            return
        end

    case 2
        Ifams  = LatticeOptData.IfamsF;
        if (length(LAT)~=length(ACHRO))
            fprintf('Warning: Incompatible input to getDVs for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nvars);
            return
        end

    case 3
        Ifams  = LatticeOptData.IfamsUC;
         if (length(LAT)~=length(UC))
            fprintf('Warning: Incompatible input to getDVs for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nvars);
            return
         end

    case 4
        Ifams  = LatticeOptData.IfamsIMC1;
         if (length(LAT)~=length(IMC1))
            fprintf('Warning: Incompatible input to getDVs for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nvars);
            return
         end    

    case 5
        Ifams  = LatticeOptData.IfamsRING;
        if (length(LAT)~=length(RING))
            fprintf('Warning: Incompatible input to getDVs for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nvars);   
            return
        end

    case 6
        Ifams  = LatticeOptData.IfamsRINGGRD;
        if (length(LAT)~=length(RINGGRD))
            fprintf('Warning: Incompatible input to getDVs for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nvars);   
            return
        end

    case 7
        Ifams  = LatticeOptData.IfamsACHROGRD;
        if (length(LAT)~=length(ACHROGRD))
            fprintf('Warning: Incompatible input to getDVs for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nvars);      
            return
        end

    otherwise
        fprintf('%s Warning: Error in getDVs, unknow lattice type nLAt = %2d \n',datetime, nLAT);
end
    
%% Reads in DVs
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

for i=1:length(bafamlist)
     Thetas = NaN(1,size(Ifams{bafamlist(i)},2));
     for l=1:size(Ifams{bafamlist(i)},1)
            Thetas(l)= LAT{Ifams{bafamlist(i)}(l)}.BendingAngle;
     end
     
     [tmax,pos] = max(abs(Thetas));
     DVs(bafamlist(i))=sign(Thetas(pos))*tmax;
end

for i=1:length(Lfamlist)
    if (not(isempty(Ifams{Lfamlist(i)})))
        DVs(Lfamlist(i))=LAT{Ifams{Lfamlist(i)}(1)}.Length;
    else
        DVs(Lfamlist(i))= NaN;
    end
end

for i=1:length(slicefamlist)
    if (not(isempty(Ifams{slicefamlist(i)})))
        DVs(slicefamlist(i))=LAT{Ifams{slicefamlist(i)}(1)}.BendingAngle;
    else
        DVs(slicefamlist(i))= NaN;
    end
end

end


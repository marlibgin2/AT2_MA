function DVs = getAllfams(nLAT, LAT, LatticeOptData)
% Gets the decision variables from a given lattice
% Uses info in strucuture LatticeOptData to define what these variables are
% Returns a one line matrix of nvars values
%
% nLAT is :
%          1 if DVs are to be taken from LAT with same structure as HACHRO in LatticeOptData
%          2 if DVs are to tbe taken from LAT with same structure as ACHRO in LatticeOptData
%          3 if DVs are to be taken from LAT with same structure as UC in LatticeOptData
%          4 if DVs are to be taken from LAT with same structure as IMC1 in LatticeOptData
%          5 if DVs are to be taken from LAT with same structure as RING in LatticeOptData
%          6 if DVS are to be taken from LAT with same structure as RINGGRD in LatticeOptData
%          7 if DVS are to be taken from LAT with same structure as ACHROGRD in LatticeOptData
%
% Note : the calling routine must check that the choice of nLAT and LAT are
% compatible with each other.

%% History
% PFT 2023,first version
% PFT 2024/08/19: included bend angles with fixed profile as set DV 

All_fams = LatticeOptData.All_fams;
nallfams = LatticeOptData.nallfams;

stdfamlist  = LatticeOptData.All_stdfamlist; 
nstdfamlist = LatticeOptData.All_nstdfamlist; 
if (isfield(LatticeOptData,'All_bafamlist'))
    bafamlist  = LatticeOptData.All_bafamlist;
else
    bafamlist = [];
end

if (isfield(LatticeOptData,'All_Lfamlist'))
    Lfamlist  = LatticeOptData.All_Lfamlist;
else
    Lfamlist = [];
end


famtype     = LatticeOptData.All_famtype;

UC           = LatticeOptData.UC;
ACHRO        = LatticeOptData.ACHRO;
HACHRO       = LatticeOptData.HACHRO;
IMC1         = LatticeOptData.IMC1;
RING         = LatticeOptData.RING;

if (isfield(LatticeOptData,'RINGGRD'))
    RINGGRD = LatticeOptData.RINGGRD;
end

if (isfield(LatticeOptData,'ACHROGRD'))
    ACHROGRD = LatticeOptData.ACHROGRD;
end

DVs = NaN(1, nallfams);

switch nLAT
    case 1
        Ifams  = LatticeOptData.IfamsAllH;
        if (length(LAT)~=length(HACHRO))
            fprintf('Warning: Incompatible input to getAllfams for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nallfams);
            return
        end
        
    case 2
        Ifams  = LatticeOptData.IfamsAllF;
        if (length(LAT)~=length(ACHRO))
            fprintf('Warning: Incompatible input to getAllfams for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nallfams);
            return
        end

    case 3
        Ifams  = LatticeOptData.IfamsAllUC;
         if (length(LAT)~=length(UC))
            fprintf('Warning: Incompatible input to getAllfams for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nallfams);
            return
         end

    case 4
        Ifams  = LatticeOptData.IfamsAllIMC1;
         if (length(LAT)~=length(IMC1))
            fprintf('Warning: Incompatible input to getAllfams for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nallfams);
            return
         end    

    case 5
        Ifams  = LatticeOptData.IfamsAllRING;
        if (length(LAT)~=length(RING))
            fprintf('Warning: Incompatible input to getAllfams for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nallfams);   
            return
        end

    case 6
        Ifams  = LatticeOptData.IfamsAllRINGGRD;
        if (length(LAT)~=length(RINGGRD))
            fprintf('Warning: Incompatible input to getAllfams for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nallfams);   
            return
        end

    case 7
        Ifams  = LatticeOptData.IfamsAllACHROGRD;
        if (length(LAT)~=length(ACHROGRD))
            fprintf('Warning: Incompatible input to getAllfams for nLAt = %2d \n',nLAT);
            DVs=NaN(1,nallfams);   
            return
        end

    otherwise
        fprintf('%s Warning: Error in getAllfams, unknow lattice type nLAt = %2d \n',datetime, nLAT);
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



end


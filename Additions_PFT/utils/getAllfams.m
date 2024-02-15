function DVs = getAllfams(nLAT, LAT, LatticeOptData)
% Gets the decision variables from a given lattice
% Uses info in strucuture LatticeOptData to define what these variables are
% Returns a one line matrix of nvars values
%
% nLAT is :
%          1 if DVs are to be taken from LAT with same structure as HACHRO in LatticeOptData
%          2 if DVs are to be taken from LAT with same structure as ACHRO in LatticeOptData
%          3 if DVs are to be taken from LAT with same structure as UC in LatticeOptData
%          4 if DVs are to be taken from LAT with same structure as IMC1 in LatticeOptData
%          5 if DVs are to be taken from LAT with same structure as RING in LatticeOptData
%
% Note : the calling routine must check that the choice of nLAT and LAT are
% compatible with each other.

All_fams = LatticeOptData.All_fams;
nallfams = LatticeOptData.nallfams;

stdfamlist  = LatticeOptData.All_stdfamlist; 
nstdfamlist = LatticeOptData.All_nstdfamlist; 
famtype     = LatticeOptData.All_famtype;

UC           = LatticeOptData.UC;
ACHRO        = LatticeOptData.ACHRO;
HACHRO       = LatticeOptData.HACHRO;
IMC1         = LatticeOptData.IMC1;
RING         = LatticeOptData.RING;

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


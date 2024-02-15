function NewLAT = setAllfams(nLAT,LAT,LatticeOptData,DVs)
%Sets decision variables to lattice LAT and returns a new lattice structure
%
% nLAT is :
%          1 if DVs are to be set to a LAT with same structure as HACHRO in LatticeOptData
%          2 if DVs are to be set to a LAT with same structure as ACHRO in LatticeOptData
%          3 if DVs are to be set to a LAT with same structure as UC in LatticeOptData
%          4 if DVs are to be taken from LAT with same structure as IMC1 in LatticeOptData
%          5 if DVs are to be taken from LAT with same structure as RING in LatticeOptData
%
% Note : the calling routine must check that the choice of nLAT and LAT are
% compatible with each other.
%   

nallfams        = LatticeOptData.nallfams;
stdfamlist      = LatticeOptData.All_stdfamlist;
nstdfamlist     = LatticeOptData.All_nstdfamlist;
famtype         = LatticeOptData.All_famtype;


UC      = LatticeOptData.UC;
ACHRO   = LatticeOptData.ACHRO;
HACHRO  = LatticeOptData.HACHRO;
IMC1    = LatticeOptData.IMC1;
RING    = LatticeOptData.RING;


switch nLAT
    case 1
        Ifams  = LatticeOptData.IfamsAllH;
         if (length(LAT)~=length(HACHRO))
            fprintf('Warning: Incompatible input to setAllfams for nLAt = %2d \n',nLAT);
            NewLAT = LAT;    
            return
        end
    case 2
        Ifams  = LatticeOptData.IfamsAllF;
        if (length(LAT)~=length(ACHRO))
            fprintf('Warning: Incompatible input to setAllfams for nLAt = %2d \n',nLAT);
            NewLAT = LAT;    
            return
        end
    case 3
        Ifams  = LatticeOptData.IfamsAllUC;
        if (length(LAT)~=length(UC))
            fprintf('Warning: Incompatible input to setAllfams for nLAt = %2d \n',nLAT);
            NewLAT = LAT;    
            return
        end
        
    case 4
        Ifams  = LatticeOptData.IfamsAllIMC1;
        if (length(LAT)~=length(IMC1))
            fprintf('Warning: Incompatible input to setAllfams for nLAt = %2d \n',nLAT);
            NewLAT = LAT;    
            return
        end
        
    case 5
        Ifams  = LatticeOptData.IfamsAllRING;
        if (length(LAT)~=length(RING))
            fprintf('Warning: Incompatible input to setAllfams for nLAt = %2d \n',nLAT);
            NewLAT = LAT;    
            return
        end
end

NewLAT = LAT;
for i=1:length(stdfamlist)
     for l=1:size(Ifams{stdfamlist(i)},1)
         if(not(isnan(DVs(stdfamlist(i)))))
            switch famtype(stdfamlist(i))
                case 1
                    NewLAT{Ifams{stdfamlist(i)}(l)}.PolynomB(1,2) = DVs(stdfamlist(i));
                    NewLAT{Ifams{stdfamlist(i)}(l)}.K = DVs(stdfamlist(i));
                case 2
                    NewLAT{Ifams{stdfamlist(i)}(l)}.PolynomB(1,3) = DVs(stdfamlist(i));
            end
         end
     end
end

for i=1:length(nstdfamlist)
    if(not(isnan(DVs(nstdfamlist(i)))))
       Ks = NaN(1,size(Ifams{nstdfamlist(i)},2));
     
       for l=1:size(Ifams{nstdfamlist(i)},1)
            Ks(l)= LAT{Ifams{nstdfamlist(i)}(l)}.PolynomB(1,2);
       end
       [kmax,pos] = max(abs(Ks));
       factor = DVs(nstdfamlist(i))/kmax*sign(Ks(pos));
     
       for l=1:size(Ifams{nstdfamlist(i)},1)
          if(not(isnan(Ks(l))))
             NewLAT{Ifams{nstdfamlist(i)}(l)}.PolynomB(1,2)= Ks(l)*factor;
             NewLAT{Ifams{nstdfamlist(i)}(l)}.K = Ks(l)*factor;
          end
       end
    end  
end
end

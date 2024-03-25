function NewLAT = setAllfamsO(nLAT,LAT,LatticeOptData,DVs)
%Sets decision variables to lattice LAT and returns a new lattice structure
%
% nLAT is :
%          1 if DVs are to be set to a LAT with same structure as HACHRO in LatticeOptData
%          2 if DVs are to be set to a LAT with same structure as ACHRO in LatticeOptData
%          3 if DVs are to be set to a LAT with same structure as UC in LatticeOptData
%          4 if DVs are to be taken from LAT with same structure as IMC1 in LatticeOptData
%          5 if DVs are to be taken from LAT with same structure as RING in LatticeOptData
%          6 if DVs are to be taken from LAT with same structure as RINGGRD in LatticeOptData
%          7 if DVs are to be taken from LAT with same structure as ACHROGRD in LatticeOptData
%
% Note : the calling routine must check that the choice of nLAT and LAT are
% compatible with each other.
%   
if (isfield(LatticeOptData,'All_fams0'))
    All_famsO = LatticeOptData.All_famsO;
    nallfamsO = LatticeOptData.nallfamsO;

    stdfamlistO  = LatticeOptData.All_stdfamlistO; 
    nstdfamlistO = LatticeOptData.All_nstdfamlistO; 
    famtypeO     = LatticeOptData.All_famtypeO;
end

nallfams        = LatticeOptData.nallfams;
stdfamlist      = LatticeOptData.All_stdfamlist;
nstdfamlist     = LatticeOptData.All_nstdfamlist;
famtype         = LatticeOptData.All_famtype;


UC      = LatticeOptData.UC;
ACHRO   = LatticeOptData.ACHRO;
HACHRO  = LatticeOptData.HACHRO;
IMC1    = LatticeOptData.IMC1;
RING    = LatticeOptData.RING;

if (isfield(LatticeOptData,'RINGGRD'))
    RINGGRD = LatticeOptData.RINGGRD;
end

if (isfield(LatticeOptData,'ACHROGRD'))
    ACHROGRD = LatticeOptData.ACHROGRD;
end

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

    case 6
        Ifams  = LatticeOptData.IfamsAllRINGGRD;
        if (length(LAT)~=length(RINGGRD))
            fprintf('Warning: Incompatible input to setAllfams for nLAt = %2d \n',nLAT);
            NewLAT = LAT;    
            return
        end

    case 7
        Ifams  = LatticeOptData.IfamsAllACHROGRDO;
        stdfamlist  = stdfamlistO; 
        nstdfamlist = nstdfamlistO;
        famtype     = famtypeO;
        if (length(LAT)~=length(ACHROGRD))
            fprintf('Warning: Incompatible input to setAllfams for nLAt = %2d \n',nLAT);
            NewLAT=LAT;   
            return
        end

    otherwise
        fprintf('%s Warning: Error in setAllfamsO, unknow lattice type nLAt = %2d \n',datetime, nLAT);

end

NewLAT = LAT;
for i=1:length(stdfamlist)
     for l=1:size(Ifams{stdfamlist(i)},1)
         if(not(isnan(DVs(stdfamlist(i)))))
            switch famtype(stdfamlist(i))
                case 1
                    NewLAT{Ifams{stdfamlist(i)}(l)}.PolynomB(1,2) = DVs(stdfamlist(i));
                    NewLAT{Ifams{stdfamlist(i)}(l)}.K = DVs(stdfamlist(i));
                case {2;3}
                    NewLAT{Ifams{stdfamlist(i)}(l)}.PolynomB(1,famtype(stdfamlist(i))+1) = DVs(stdfamlist(i));
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

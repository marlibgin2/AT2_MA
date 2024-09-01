function NewLAT = setAllfams(nLAT,LAT,LatticeOptData,DVs)
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
   
%% History
% PFT 2023,first version
% PFT 2024/08/19: included bend angles with fixed profile as DV 
% PFT 2024/08/24: included distance as DV
% PFT 2024/08/28: included handling of element slice bend angles as DV
%
%%
nallfams        = LatticeOptData.nallfams;
stdfamlist      = LatticeOptData.All_stdfamlist;
nstdfamlist     = LatticeOptData.All_nstdfamlist;
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

if (isfield(LatticeOptData,'All_slicefamlist'))
    slicefamlist = LatticeOptData.All_slicefamlist;
else
    slicefamlist = [];
end

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
        Ifams  = LatticeOptData.IfamsAllACHROGRD;
        if (length(LAT)~=length(ACHROGRD))
            fprintf('Warning: Incompatible input to setAllfams for nLAt = %2d \n',nLAT);
            NewLAT = LAT;    
            return
        end

    otherwise
        fprintf('%s Warning: Error in setAllfams, unknow lattice type nLAt = %2d \n',datetime, nLAT);
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

for i=1:length(bafamlist)
    if(not(isnan(DVs(bafamlist(i)))))
       Thetas = NaN(1,size(Ifams{bafamlist(i)},2));
     
       for l=1:size(Ifams{bafamlist(i)},1)
            Thetas(l)= LAT{Ifams{bafamlist(i)}(l)}.BendingAngle;
       end
       [tmax,pos] = max(abs(Thetas));
       factor = DVs(bafamlist(i))/tmax*sign(Thetas(pos));
     
       for l=1:size(Ifams{bafamlist(i)},1)
          if(not(isnan(Thetas(l))))
             NewLAT{Ifams{bafamlist(i)}(l)}.BendingAngle = Thetas(l)*factor;
          end
       end
    end  
end

for i=1:length(Lfamlist)
     for l=1:size(Ifams{Lfamlist(i)},1)
         if(not(isnan(DVs(Lfamlist(i)))))
             NewLAT{Ifams{Lfamlist(i)}(l)}.Length = DVs(Lfamlist(i));
         end
     end
end

for i=1:length(slicefamlist)
     for l=1:size(Ifams{slicefamlist(i)},1)
         if(not(isnan(DVs(slicefamlist(i)))))
             NewLAT{Ifams{slicefamlist(i)}(l)}.BendingAngle = DVs(slicefamlist(i));
         end
     end
end

end

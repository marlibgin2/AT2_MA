function NewLAT = setLins(nLAT,LAT,LatticeOptData,DVs)
%Sets decision variables of quadrupole gradients to lattice LAT 
% and returns a new lattice structure
%
% nLAT is :
%          1 if DVs are to be set to a LAT with same structure as HACHRO in LatticeOptData
%          2 if DVs are to be set to a LAT with same structure as ACHRO in LatticeOptData
%          3 if DVs are to be set to a LAT with same structure as UC in LatticeOptData
%          4 if DVs are to be taken from LAT with same structure as IMC1 in LatticeOptData
%          5 if DVs are to be taken from LAT with same structure as RING in LatticeOptData
%          6 if DVS are to be taken from LAT with same structure as RINGGRD in LatticeOptData
%          7 if DVS are to be taken from LAT with same structure as ACHROGRD in LatticeOptData
%
% Note : the calling routine must check that the choice of nLAT and LAT are
% compatible with each other.
%   
optMode = LatticeOptData.optMode;
%
% check for backward compatibility
%
%
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

UC           = LatticeOptData.UC;
ACHRO        = LatticeOptData.ACHRO;
HACHRO       = LatticeOptData.HACHRO;
RING         = LatticeOptData.RING;

if (isfield(LatticeOptData,'IMC1'))
    IMC1 = LatticeOptData.IMC1;
end

if (isfield(LatticeOptData,'RINGGRD'))
    RINGGRD = LatticeOptData.RINGGRD;
end

if (isfield(LatticeOptData,'ACHROGRD'))
    ACHROGRD = LatticeOptData.ACHROGRD;
end

switch nLAT
    case 1
        Ifams  = LatticeOptData.IfamsHlin;
         if (length(LAT)~=length(HACHRO))
            fprintf('Warning: Incompatible input to setLins for nLAt = %2d \n',nLAT);
            NewLAT = LAT;    
            return
         end

    case 2
        Ifams  = LatticeOptData.IfamsFlin;
        if (length(LAT)~=length(ACHRO))
            fprintf('Warning: Incompatible input to setLins for nLAt = %2d \n',nLAT);
            NewLAT = LAT;    
            return
        end

    case 3
        Ifams  = LatticeOptData.IfamsUClin;
        if (length(LAT)~=length(UC))
            fprintf('Warning: Incompatible input to setLins for nLAt = %2d \n',nLAT);
            NewLAT = LAT;    
            return
        end

        
    case 4
        Ifams  = LatticeOptData.IfamsIMC1lin;
        if (length(LAT)~=length(IMC1))
            fprintf('Warning: Incompatible input to setLins for nLAt = %2d \n',nLAT);
            NewLAT = LAT;    
            return
        end
        
     case 5
        Ifams  = LatticeOptData.IfamsRINGlin;
        if (length(LAT)~=length(RING))
            fprintf('Warning: Incompatible input to setLins for nLAt = %2d \n',nLAT);
            NewLAT = LAT;    
            return
        end

     case 6
        Ifams  = LatticeOptData.IfamsRINGGRDlin;
        if (length(LAT)~=length(RINGGRD))
            fprintf('Warning: Incompatible input to setLins for nLAt = %2d \n',nLAT);
            NewLAT = LAT;    
            return
        end

     case 7
        Ifams  = LatticeOptData.IfamsACHROGRDlin;
        if (length(LAT)~=length(ACHROGRD))
            fprintf('Warning: Incompatible input to getLins for nLAt = %2d \n',nLAT);
            NewLAT = LAT;     
            return
        end

    otherwise
        fprintf('%s Warning: Error in setLins, unknow lattice type nLAt = %2d \n',datetime, nLAT);
end

NewLAT = LAT;
for i=1:length(stdfamlist)
     for l=1:size(Ifams{stdfamlist(i)},1)
         if(not(isnan(DVs(stdfamlist(i)))))
            switch famtype(stdfamlist(i))
                case 1
                    NewLAT{Ifams{stdfamlist(i)}(l)}.PolynomB(1,2) = DVs(stdfamlist(i));
                    NewLAT{Ifams{stdfamlist(i)}(l)}.K = DVs(stdfamlist(i));
                case {2,3}
                    NewLAT{Ifams{stdfamlist(i)}(l)}.PolynomB(1,famtype(stdfamlist(i))) = DVs(stdfamlist(i));
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


function OBJ = nOBJfcn_m4U(X, rin)
verbo = true;

% -----------------------------------------
% alter the lattice according to parameters
% -----------------------------------------


rrr  = alter_m4U_lattice(X, rin, 'B');
s = findspos(rrr,1:length(rrr)+1); se = s(end);
% guess the periodicity
P = round(528/se);      % <MSj comment> I've got a function 'estimateperiodicity' that does this via FFT. Likely can use a persistent variable re-initializing after 5 seconds of no function call as it's not going to change.

% ----------------------------
% force chromaticity to be 1/1
% ----------------------------
changechro=1;
if changechro == 1
    try
        if verbo; disp('OBJ -- fitting chromaticity to [1,1]'); end
        rrr = atfitchrom(rrr,[1/P 1/P],'S1','S2');
    catch
        if verbo; disp('cannot fit chromaticity ... machine unstable'); end
    end
end

object = 'emix_DA';
object = 'DA_TuneSpread';
switch object
    case 'emix_DA'

        try
            ATS = atsummary_ring(rrr);
            if isnan(ATS.naturalEmittance)
                error('no emittance -- retry')
            end
        catch
            ATS.naturalEmittance = 1e19;
        end
        emix = ATS.naturalEmittance;

        if isnan(emix)
            emix = 1e19;
        end

        r0     = 3e-3;
        nsteps = 20; %61;
        nturns = 250;
        dp     = 0.0;
        res    = 0.25e-3;
        alpha  = 1.1;
        DA     = modelDA_sim_par(rrr, r0, nsteps, nturns*P, dp, res, alpha); %%with parfor
        % DA     = modelDA_sim(rrr, r0, nsteps, nturns*P, dp, res, alpha);     %% with for

        %
        % clip tall solutions
        % 05012024
        %
        % for i=1:nsteps
        %     if DA(i,2)>5e-3
        %         DA(i,2) = 5e-3;
        %     end
        % end

        DA_Area = -polyarea(DA(:,1),DA(:,2));

        if isnan(DA_Area)
            DA_Area = 1e19;
        end

        % try
        %     DA_Area = -calcDA_Area(DA);
        % catch
        %     DA_Area = 1e19;
        % end

        disp(['emix = ' num2str(emix*1e12,8) ' pmrad, DA_Area = ' num2str(-DA_Area*1e6,8) ' mm²'])
        disp(['X = [' num2str(X,10) ']'])
        writelines(['emix = ' num2str(emix*1e12,8) ' pmrad, DA_Area = ' num2str(-DA_Area*1e6,8) ' mm²'],"temp.txt",WriteMode="append")
        writelines(['X = [' num2str(X,10) ']'],"temp.txt",WriteMode="append")
        writelines(['                          '],"temp.txt",WriteMode="append")
        writelines(['                          '],"temp.txt",WriteMode="append")

        OBJ = [emix*1e12, DA_Area*1e6];     % Rescaling to pm rad and mm^2

    case 'DA_TuneSpread'
        r0     = 3e-3;
        nsteps = 20; %61;
        nturns = 250;
        dp     = 0.0;
        res    = 0.25e-3;
        alpha  = 1.1;
        DA     = modelDA_sim_par(rrr, r0, nsteps, nturns*P, dp, res, alpha); %%with parfor
        DA_Area = -polyarea(DA(:,1),DA(:,2));

        if isnan(DA_Area)
            DA_Area = 1e19;
        end

%
% [x0pos, y0pos, nuxpos, nuypos, diffuvec, WA] = DA_NU_FMA_par(rrr,'fma',2000*P,12,8,12e-3,7e-3,'newcalc','aurora cluster','cluster cores',32);
        [x0pos, y0pos, nuxpos, nuypos, diffuvec, WA] = DA_NU_FMA_par(rrr,'fma',2000,24,16,12e-3,7e-3,'newcalc','cluster cores',32);
        x0 = mean(x0pos,1); y0 = mean(y0pos,1);
        x0=x0*1e3; y0=y0*1e3;
        sele1 = find(abs(x0)<15&abs(y0)<15);
        %x0=x0(sele1); y0=y0(sele1); 
        DAmm = DA*1e3;
        [in,on] = inpolygon(x0,y0,DAmm(:,1),DAmm(:,2));
        sele2 = find(in>0);
        sele = intersect(sele1, sele2);
        nux = mean(nuxpos,1); nuy = mean(nuypos,1);
        nux = nux(sele); nuy = nuy(sele);
        % nux=nux(nux>0&&nux<0.5&&nuy>0&&nuy<0.5) 
        % avoid points coalescing to (0,0)
        snux = std(nux); snuy=std(nuy); snu=sqrt(snux^2+snuy^2);
        if isnan(snu)
            snu = 1e19;
        end

        OBJ = [DA_Area*1e6, snu];     % Rescaling to pm rad and mm^2

end
end

function DA_Area = calcDA_Area(DA)
a = []; 
A = [0, 0];
for i = 1:length(DA)-1
    B = DA(i,  :);
    C = DA(i+1,:);
    a(i) = (A(1) * (B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2)))/2;
end
%W = 1 - abs(DA(1,1) + DA(end,1))/(DA(1,1) - DA(end,1)); % weight to balance the X-asymmetry observed in some cases 
%Lx = (DA(1,1)-DA(end,1))/2; % average width (case of nstep=9 == 10 pts)
%Ly = (DA(5,2)+DA(6,2))/2; % average height
%Lx = (DA(1,1)-DA(end,1))/2; % average width (case of nstep=14 == 15pts
%Ly = (DA(8,2)); % average height
%AspRa = abs(Lx/Ly); % Aspect Ratio (want a DA large at the base, not too pointy)
DA_Area = sum(a);% * W;% * AspRa;
end


function ScanResults = glass_scan(Ns, Ntmax, Ni, MaxEmit, MaxChro, MaxEta, savef, choosef, folder)
%GLASS_SCAN Performs GLASS scan of 3 GeV ring 7 BA
%   inputs:
%          Ns : number of steps for each one of the 7 scan variables. Total # points = Ns^7
%          Ntmax : maximum number of scans to perform. Allows stopping the loop at 
%                  a point before the total # of scans is realized and save
%                  files in intermediate files
%          Ni : index of first scan - allows resuming an interrupted large scan.
%          MaxEmit : only lattices with emittances larger than zero ans
%                    smaller than MAX Emit are recorded.
%          MaxChro : maximum absolute value of both chromaticities.
%          MaxEta  : maximum absolute value of dispersion at centre of long
%                    sraight.
%          savef   : if 'Y' , data is saved on a file in the form of a
%                    structure that is also returned by the function
%          choosef : if 'N', all stable (abs(Tr M)<2) lattices are recorded. If 'Y', only
%                    lattices with 0<Jx<3, 0<Emit<MaxEmit, abs(ChrosX,y)< MAx Chro,
%                    abs(Eta at centre of long straight) < Max Eta.
%          folder:   folder in whihc files will be saved in any. The folder
%                    must exist and be in the matlab path before the function is called. When using
%                    the LongScan function to call glass_scan, the folder
%                    is input through a dialog.
%
% Outputs: Structure (possiby sved onto a file) with the following fields
%
%         scan_fams: cell array of srrings - magnet families used in the
%                    scan
%         T_rb_deg  : bending angle of reverse bend in the unit cell (deg), 
%         
%         Ns: number of steps for each one of the 7 scan variables.
%         Maxs:maximum of each variable (m**-2)          
%         Mins: minium of each variable (m**-2)
%         deltas: step ineach variable  (m**-´2)
%         Ni: index of first scan
%         MaxEmit: maximum emittance (pmrad)
%         MaxChro: maximum cromaticities
%         MaxEta: maximum dispersion at centre of long stright (m)
%         choosef: 'N': does not select stable latuces to record
%         nstab: number of stable lattices
%         nrecord: number of recorded lattices
%         nerror : number of lattices that generated erors in atlinopt4
%         Values : Matrix collecting all recorder lattices. One line per lattice wth the folowing columns
%           1: qfend (m**-2)
%           2: qdend (m**-2)
%           3: reversebend (m**-2)
%           4: diph (m**-2)
%           5: dipm (m**-2)
%           6: reversebendmc1 (m**-2)
%           7: reversebendmc2 (m**-2)
%           8: Emittance (pmrad)
%           9: Jx
%          10: TuneX (fractional part one achromat)
%          11: TuneY (fractional part one achromat)
%          12: ChroX
%          13: ChroY
%          14: Etax (m)
%          15: Qx_ring
%          16: Qy_ring
%
% The global variables
% below must be initlized before calling the function - this can be done with the function 
% max4_UpgradeStudies_20220305_AT2
%
% This function is best used from he function LongScan that automates
% caling it repeatedily so data can be saved in intermediate steps.
%

global HACHRORBSSIM % half achromat supersimpified (only one slice per dipole, concatenatd drifts) 
global ACHRORBSSIM  % complete simplified achromat


nscan = 0;
nfrac = 0;
nstab = 0;
nrecord = 0;
nerror = 0;

NE=length(HACHRORBSSIM); % 
props=atCheckRingProperties(HACHRORBSSIM);

% Documents the reverse bend angle
%
I_rb = find(atgetcells(HACHRORBSSIM, 'FamName','reversebend_sim'));
T_rb = atgetfieldvalues(HACHRORBSSIM, I_rb(1), 'BendingAngle');
T_rb = T_rb*180/pi;

% Magnet families used in the scan
%
scan_fams = {'qfend_sim';...
             'qdend_sim';
             'reversebend_sim';...
             'sshdip';...
             'ssdipm';...
             'reversebendmc1_sim';...
             'reversebendmc2_sim'};
Nsteps   = [1;1;1;1;1;1;1]*Ns;         
Maxs     = [5.00;0.0;5.00;-0.80;-0.80; 5.00; 5.00];
Mins     = [0.00;-5.0;0.00;-1.25;-1.25; 0.00; 0.00];
deltas   = (Maxs-Mins)./(Nsteps-1);

nvars = size(scan_fams,1);

nt=Nsteps(1);

for i=2:nvars
    nt = nt*Nsteps(i);
end

for i=1:nvars
    I_famsH{i} =  find(atgetcells(HACHRORBSSIM,'FamName',scan_fams{i}));
end

for i=1:nvars
    I_famsF{i} =  find(atgetcells(ACHRORBSSIM,'FamName',scan_fams{i}));
end

angle=atgetfieldvalues(ACHRORBSSIM,'BendingAngle');
isdipole=isfinite(angle) & (angle~=0);

fprintf('Total number of evaluations  : %12d  \n', nt);
nt=min(nt-Ni+1,Ntmax);
fprintf('Actual number of evaluations : %12d  \n', nt);

Vals = zeros(nt,nvars+9);
fb=waitbar(0,'Starting Calculation...');

LAT=HACHRORBSSIM;
t_est = nt/10000*5;

fprintf('Estimated time : %5.2f seconds \n', t_est);
Kval=zeros(1,nvars);
tic;
for i=Ni:Ni+nt-1
    for j=1:nvars
        k = mod(floor((i-1)/(Nsteps(j)^(nvars-j))),Nsteps(j));
        Kval(j) = Mins(j)+k*deltas(j);
        for l=1:size(I_famsH{j},1)
            LAT{I_famsH{j}(l)}.PolynomB(1,2) = Kval(j);
            LAT{I_famsH{j}(l)}.K = Kval(j);
        end
    end
    twissdatah = findm44_fast(LAT,NE,props);
    TRx = 2*(twissdatah(1,1)*twissdatah(2,2)+twissdatah(1,2)*twissdatah(2,1));
    TRy = 2*(twissdatah(3,3)*twissdatah(4,4)+twissdatah(3,4)*twissdatah(4,3));
    
    if ((abs(TRx)<1.98)&&(abs(TRy)<1.98))
         nstab=nstab+1;
         LAT=ACHRORBSSIM;
         for j=1:nvars
             for l=1:size(I_famsF{j},1)
                 LAT{I_famsF{j}(l)}.PolynomB(1,2) = Kval(j);
                 LAT{I_famsF{j}(l)}.K = Kval(j);
             end
         end
         try
            rpara = atsummary_fast(LAT,isdipole);
            Emitt = rpara.naturalEmittance*1E12;
            Jx    = rpara.damping(1);
            ChroX = rpara.chromaticity(1);
            ChroY = rpara.chromaticity(2);
            Etax  = rpara.etax;
            tuneX = rpara.tunes(1);
            tuneY = rpara.tunes(2);
            Qx_ring = rpara.Qx_ring;
            Qy_ring = rpara.Qy_ring;
            
            if ( (strcmp(choosef,'N'))||((Emitt>0) && (Emitt<MaxEmit) && (Jx<3) && ...
               abs(ChroX)<MaxChro && abs(ChroY)<MaxChro && abs(Etax)<MaxEta))
              nrecord=nrecord+1;
              Vals(nrecord,1:nvars)= Kval;
              Vals(nrecord,nvars+1)= Emitt;
              Vals(nrecord,nvars+2)= Jx;
              Vals(nrecord,nvars+3)= tuneX;
              Vals(nrecord,nvars+4)= tuneY;
              Vals(nrecord,nvars+5)= ChroX;
              Vals(nrecord,nvars+6)= ChroY;
              Vals(nrecord,nvars+7)= Etax;
              Vals(nrecord,nvars+8)= Qx_ring;
              Vals(nrecord,nvars+9)= Qy_ring;
            end
         catch ME
            nerror=nerror+1;
            nrecord=nrecord+1;
            fprintf('%s \n', '***************************************');
            fprintf('Error executing atsummary at i= %d11 \n', i);
            fprintf('Error identifier = %s \n', ME.identifier);
            fprintf('Traces = %6.2f %6.2f \n', TRx, TRy);
  %          for m=1:nvars
  %              fprintf(strcat(scan_fams{m} ,' = %8.4f \n'), Kval(m));
  %         end
            
            Vals(nrecord,1:nvars)= Kval;
            Vals(nrecord,nvars+1)= NaN;
            Vals(nrecord,nvars+2)= NaN;
            Vals(nrecord,nvars+3)= NaN;
            Vals(nrecord,nvars+4)= NaN;
            Vals(nrecord,nvars+5)= NaN;
            Vals(nrecord,nvars+6)= NaN;
         end
    end
    nscan=nscan+1;
    nfracnew=nscan/nt*100;
    if ((nfracnew-nfrac)>10)
        waitbar(nscan/nt,fb,...
                 strcat(sprintf('%3.0f %s',nfracnew,'%'),...
                        sprintf('%6d stable lat.', nstab),...
                        sprintf('%4d recorded', nrecord),...
                        sprintf('%3d errors', nerror)));
                    
        nfrac=nfracnew;
    end
end
delete(fb);
fprintf('%7d stable lattices found  \n', nstab);
fprintf('%5d lattices with errors  \n', nerror);
fprintf('%5d lattices recorded \n', nrecord);
%fprintf('Emit <  %4.2f  pmrad \n', MaxEmit);
%fprintf('|Chroms| <  %4.2f  \n', MaxChro);
%fprintf('|Eta x| <  %4.2f  mm \n', MaxEta*1000);
%fprintf('Jx < %1d \n', 3);

filename=strcat(folder,'\GScan_',datestr(now,30));
Vals(nrecord+1:end,:)=[];
ScanResults=struct('scan_fams',{scan_fams},'T_rb_deg', T_rb,'nscan',nscan,...
                   'Ns', {Ns},'Maxs', Maxs, 'Mins', Mins, 'deltas',deltas,...
                   'Ni', {Ni},'nstab',{nstab},'nrecord', nrecord,...
                   'nerror', nerror, 'MaxEmit', MaxEmit, ...
                   'MaxChro', MaxChro, 'MaxEta', MaxEta,...
                   'choosef', choosef,...
                   'Values',Vals);
               
if (strcmp(savef,'Y'))
       save(filename,'ScanResults');
end
toc;
end


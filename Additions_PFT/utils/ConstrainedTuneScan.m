function CTScan = ConstrainedTuneScan(ACHRO, LatticeOptData, X0Range, Qrange, qrange, Npq, plottunescanf, verbose)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    chrom_fams=LatticeOptData.chrom_fams;
    isdipole = LatticeOptData.isdipole;
    DAoptions = LatticeOptData.DAoptions;
    chroms0 = DAoptions.chroms0;
    Nitchro = DAoptions.Nitchro;
    TolChrom = DAoptions.TolChrom;
    ACHRO_scan=ACHRO;
    Qxmin = Qrange(1);
    Qxmax = Qrange(2);
    Qymin = Qrange(3);
    Qymax = Qrange(4);

    qxmin = qrange(1);
    qxmax = qrange(2);
    qymin = qrange(3);
    qymax = qrange(4);

    Npqx = Npq(1);
    Npqy = Npq(2);

    nx = (Qxmax-Qxmin+1)*Npqx;
    ny = (Qymax-Qymin+1)*Npqy;
    dqx = (qxmax-qxmin)/(Npqx-1);
    dqy = (qymax-qymin)/(Npqy-1);

    ntot = nx*ny;
    fprintf('Starting tune scan with %3d grid points \n', ntot);
    %fb=waitbar(0,'Starting Tune Scan...', 'Name','Tune Scan Progress', 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    %setappdata(fb,'canceling',0);
    I_sc1 = find(atgetcells(ACHRO, 'FamName', chrom_fams{1}));
    I_sc2 = find(atgetcells(ACHRO, 'FamName', chrom_fams{2}));

    DAScan = nan(nx,ny);
    EmitScan = nan(nx,ny);
    BetaX0Scan = nan(nx,ny);
    BetaY0Scan = nan(nx,ny);
    JxScan     = nan(nx,ny);
    qxscan     = nan(1,nx);
    qyscan     = nan(1,ny);

    
    k=1;
    for i=Qxmin:Qxmax
        for j=1:Npqx
            qxscan(k)=i+qxmin+dqx*(j-1);
            k=k+1;
        end
    end

    k=1;
    for i=Qymin:Qymax
        for j=1:Npqy
            qyscan(k)=i+qymin+dqy*(j-1);
            k=k+1;
        end
    end

    
    nfrac=0;
    nscan=0;

    for i=1:nx
        Qxfit = qxscan(i);
        for j=1:ny
            Qyfit = qyscan(j);
            try
                [ACHRO_scan, Penalty]= LinMatch(ACHRO_scan,'Y','Y',[Qxfit, Qyfit],'N',[0 100],30,[],[],500, LatticeOptData,'N');
                if (strcmpi(verbose,'Y'))
                    fprintf('Ring Tunes fit to [ %4.2f , %4.2f ] with penalty = %6.2e \n', Qxfit, Qyfit, Penalty);
                end
            catch ME
                fprintf('ExMOGA Error in tune fit at ring tunes  [ %4.2f , %4.2f ] \n', Qxfit, Qyfit);
                fprintf('Error message was:%s \n',ME.message);
                ACHRO_scan=ACHRO;
                EmitScan(i,j)   = NaN;
                BetaX0Scan(i,j) = NaN;
                BetaY0Scan(i,j) = NaN;
                JxScan(i,j)     = NaN;
                DAScan(i,j)     = NaN;
                nscan=nscan+1;
                nfracnew=nscan/ntot*100;
                if ((nfracnew-nfrac)>1)
                    waitbar(nscan/ntot,fb,strcat(sprintf('%3.0f %s',nfracnew,'%')));
                    nfrac=nfracnew;
                end
                continue
            end
            try
                rparascan=atsummary_fast(ACHRO_scan,isdipole);
                EmitScan(i,j)   = rparascan.naturalEmittance;
                BetaX0Scan(i,j) = rparascan.beta0(1);
                BetaY0Scan(i,j) = rparascan.beta0(2);
                JxScan(i,j)     = rparascan.damping(1);
            catch ME
                fprintf('ExMOGA Error in atsummary at ring tunes  [ %4.2f , %4.2f ] \n', Qxfit, Qyfit);
                fprintf('Error message was:%s \n',ME.message);
                EmitScan(i,j)   = NaN;
                BetaX0Scan(i,j) = NaN;
                BetaY0Scan(i,j) = NaN;
                JxScan(i,j)     = NaN;
                DAScan(i,j)     = NaN;
                nscan=nscan+1;
                nfracnew=nscan/ntot*100;
                if ((nfracnew-nfrac)>1)
                    waitbar(nscan/ntot,fb,strcat(sprintf('%3.0f %s',nfracnew,'%')));
                    nfrac=nfracnew;
                end
                continue
            end
            %
            % Fits chromaticity 
            %
            try 
                [ACHRO_scan, Penalty, its]=fitchroit(ACHRO_scan, chrom_fams, chroms0, Nitchro, TolChrom); 
                K_sc1  = atgetfieldvalues(ACHRO_scan, I_sc1, 'PolynomB', {3});
                K_sc2  = atgetfieldvalues(ACHRO_scan, I_sc2, 'PolynomB', {3});
                Sc1    = K_sc1(1);
                Sc2    = K_sc2(1);
            if (strcmp(verbose,'Y'))
                fprintf('Chromaticity matched with penalty = %6.2e in %2d iterations\n', Penalty, its);
            end
            catch ME
                fprintf('Error in ExMOGA: chromaticity fit \n');
                fprintf('Error message was: %s \n',ME.message);
                LAT_scan = LAT_tune;
                Sc1=NaN;
                Sc2=NaN;
                DAScan(i,j) = NaN;
                nscan=nscan+1;
                nfracnew=nscan/ntot*100;
                if ((nfracnew-nfrac)>1)
                    waitbar(nscan/ntot,fb,strcat(sprintf('%3.0f %s',nfracnew,'%')));
                    nfrac=nfracnew;
                end
                continue
            end
              
            if (not(isnan(Sc1)&&not(isnan(Sc2))))
                DA=CalcPlotDA(LAT_scan,DAoptions,'N');
            else
                DA=NaN;
            end
            DAScan(i,j)=DA;
            if (strcmp(verbose,'Y'))
                fprintf('i = %3d j = %3d qx= %5.2f qy= %5.2f DA = %5.2f mm**2 Emit = %5.2f pmrad \n', i, j, Qxfit, Qyfit, DA, EmitScan(i,j)*1e12);
            end
            nscan=nscan+1;
            nfracnew=nscan/ntot*100;
            if ((nfracnew-nfrac)>1)
                    waitbar(nscan/ntot,fb,strcat(sprintf('%3.0f %s',nfracnew,'%')));
                    nfrac=nfracnew;
            end
            if getappdata(fb,'canceling')
                break
            end
        end
        if getappdata(fb,'canceling')
           fprintf ('Canceling at nscan = %5d i = %3d j = %3d qx= %5.2f qy= %5.2f \n', nscan, i, j, Qxfit, Qyfit);
           break
        end
    end 
    delete(fb);
    CTScan.qxscan=qxscan;
    CTScan.qyscan=qyscan;
    CTScan.DAScan=DAScan;
    CTScan.EmitScan=EmitScan*1e12;
    CTScan.BetaX0Scan = BetaX0Scan;
    CTScan.BetaY0Scan = BetaY0Scan;
    CTScan.JxScan = JxScan;
    if (strcmp(plottunescanf,'Y'))
        PlotEMTuneScan(CTScan);
    end
    [mdaj,jmaxda] = max(DAScan');[maxda,imaxda] = max(mdaj);
    CTScan.DAmax = maxda;
    CTScan.qxmaxDA = qxscan(imaxda);
    CTScan.qymaxDA = qyscan(jmaxda(imaxda));
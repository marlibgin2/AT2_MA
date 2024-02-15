function res = RBTuneScan_AT2(LatticeOptData,nuc,fixedThetaf, dTheta, fitchromf, verbosef,qx0, qy0, dqx, dqy, Nqx, Nqy )
%RBTuneScan: Scans unit cell tunes and calculates beam parameters for
%   a ring composed only of nuc unit cells. If chosen, chromaticity is fit to 1/120 for
%   both planes, corresponding to +1 chromaticity for a full ring.
%  
if (strcmp(verbosef,'Y'))
    fprintf('----- \n');
    fprintf('%s Starting UC Tune Scan \n', datetime);
    fprintf('----- \n');
end
tic;
if (strcmp(fixedThetaf,'Y'))
    dTheta     = LatticeOptData.Trb;
end

res.dTheta = dTheta;
res.nuc    = nuc;
res.qx0    = qx0;
res.qy0    = qy0;
res.dqx    = dqx;
res.dqy    = dqy;
res.Nqx    = Nqx;
res.Nqy    = Nqy;

res.data=[];
formatSpec = 'i = %3d  j = %3d  qx = %8.3f  qy = %12.7f  \n';
for i=-Nqx:Nqx
    qx = qx0+i*dqx;
    for j=-Nqy:Nqy
        qy = qy0 + j*dqy;                 
        tunes = [qx qy];
        if (strcmp(fitchromf,'Y'))
            [pars,UC] = RBTuneParam_AT2(nuc, tunes, LatticeOptData, fixedThetaf, dTheta, 'N', [1/120 1/120]);
        else
            [pars,UC] = RBTuneParam_AT2(nuc, tunes, LatticeOptData, fixedThetaf, dTheta, 'N');
        end
        if (pars(12)>=3.0) % Checks for longitudinally unstable motion (Jx>3)
            pars(3:56)=NaN;
        end
        LatticeOptData.UC = UC;
        res.data=cat(1,res.data,real(pars));
        if (strcmp(verbosef,'Y'))
            fprintf(formatSpec,i,j,qx,qy);
        end
    end
end
toc;
end



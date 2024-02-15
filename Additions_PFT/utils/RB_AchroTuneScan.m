function TuneScanDet = RB_AchroTuneScan(LatticeOptData,UCTuneTable,...
                    verbosef,verboselevel,nit,TolTunes,TolDisp,...
                    TolAlphas, TolX, lb, ub)
%
%Performs matching of dispersion and alphas for 
%a list of desired unit cell tunes

ntunes = size(UCTuneTable,1);

TuneScanDet = cell(ntunes,1);
for i=1:ntunes
    uctunex=UCTuneTable(i,1);
    uctuney=UCTuneTable(i,2);
    if(not(isnan(uctuney)&&not(isnan(uctunex))))
        rp = RBAchroTuneParam_AT2(LatticeOptData, 'Y', ...
                 0.0, verbosef, verboselevel, 'N', 'N', ...
                 'Y', [UCTuneTable(i,2) UCTuneTable(i,1)], nit, ...
                 TolTunes, TolDisp, TolAlphas, TolX, lb, ub);
        TuneScanDet{i}=rp; 
    else
        TuneScanDet{i} = NaN;
    end
end


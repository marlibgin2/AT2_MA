function uctunex = uctunexvsemit(LatticeOptData,uctuney,emittance,nit,TolTunes,verbosef)
%finds the required horizontal unit cell phase advance to reach a desired 
% emittance at a given vertical unit cell phase advance
%
UC=LatticeOptData.UC;
%uctunex=emitvsuctunes(LatticeOptData,[0.05 uctuney],nit,TolTunes,emittance)
%uctunex=emitvsuctunes(LatticeOptData,[0.45 uctuney],nit,TolTunes,emittance)
%uctunex=emitvsuctunes(LatticeOptData,[3/7 1/7],nit,TolTunes,emittance)
if (strcmp(verbosef,'Y'))
    options = optimset('Display','iter','TolX',1E-5);
else
    options = optimset('Display','off','TolX',1E-5);
end


try
    uctunex = fzero(@(x)emitvsuctunes...
               (LatticeOptData,[x uctuney],nit,TolTunes,emittance),...
               0.25,options);
catch ME
    fprintf('Problems finding horizontal tune \n');
    fprintf ('Error = %s \n', ME.message); 
    uctunex=NaN;
end
end

function emit=emitvsuctunes(LatticeOptData, uctunes, nit, TolTunes,emittance)

UC=LatticeOptData.UC;
isdipoleUC = LatticeOptData.isdipoleUC;
scan_fams=LatticeOptData.scan_fams;
try
  [UC_T, its, penalty, ftunes]=fittuneRS(UC,uctunes,scan_fams{3},...
                               scan_fams{4},nit, TolTunes,'No');
   rp=atsummary_fast(UC_T,isdipoleUC); 
%   rp=atsummary(UC_T); 
   
   emit=rp.naturalEmittance-emittance;
%   emit=rp.naturalEmittance;
catch
    fprintf('Problems fitting tunes %8.3f %8.3f \n', uctunes(1),uctunes(2));
    emit=NaN;
end


end


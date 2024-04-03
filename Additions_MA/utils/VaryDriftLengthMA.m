function R=VaryDriftLengthMA(R,DL,fam1,fam2,SumD)
%
% changes length of two drifts families, keeping total length constant
% initially thoght to move QM1 into the long straight 
%
id1 = findcells(R,'FamName',fam1); d1 = getcellstruct(R,'Length',id1); d1=2*d1(1);
id2 = findcells(R,'FamName',fam2); d2 = getcellstruct(R,'Length',id2); d2=d2(1);
disp(['SD_before =' num2str(d1+d2) ' (m)'])
%nd1=length(id1);
%nd2=length(id2);
R=setcellstruct(R,'Length',id1,DL*ones(size(id1)));
R=setcellstruct(R,'Length',id2,(SumD-2*DL)*ones(size(id2)));
d1 = getcellstruct(R,'Length',id1); d1=2*d1(1);
d2 = getcellstruct(R,'Length',id2); d2=d2(1);
disp(['SD_after =' num2str(d1+d2) ' (m)'])
return
function QFMultErrorScan(ACHRO,an,bn,r0,n0,famnames,DAoptions)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
AC = ACHRO;
CalcPlotDA(AC,DAoptions,'plot');
hold;
AC=AddMultipoleError(ACHRO,an,bn,r0,famnames,n0,1);
[DA,DAV]=calcDA_fast(AC,DAoptions,0.0);
plot(DAV(:,1)*1000,DAV(:,2)*1000,'r-o');
AC=AddMultipoleError(ACHRO,an,bn,r0,famnames,n0,10);
[DA,DAV]=calcDA_fast(AC,DAoptions,0.0);
plot(DAV(:,1)*1000,DAV(:,2)*1000,'g-o');
AC=AddMultipoleError(ACHRO,an,bn,r0,famnames,n0,100);
[DA,DAV]=calcDA_fast(AC,DAoptions,0.0);
plot(DAV(:,1)*1000,DAV(:,2)*1000,'k-o');
AC=AddMultipoleError(ACHRO,an,bn,r0,famnames,n0,1000);
[DA,DAV]=calcDA_fast(AC,DAoptions,0.0);
plot(DAV(:,1)*1000,DAV(:,2)*1000,'c-o');
legend('0','1','10','100','1000');
end
function [s1,s2,s3,s4]=PlotEMTuneScan(rm)
% Plots resuts of two-dimensional tune scan - data generated, e.g. from ExMOGA
% Matrices are given in the inpout structure rm
%
qxscan = rm.qxscan;
qyscan = rm.qyscan;
DAScan = rm.DAScan;
BetaX0Scan = rm.BetaX0Scan;
BetaY0Scan = rm.BetaY0Scan;
EmitScan = rm.EmitScan;

figure;s1=pcolor(qxscan,qyscan,DAScan');s1.FaceColor='interp';s1.EdgeColor='none';xlabel('Qx');ylabel('Qy');title('DA [mm**2]');colorbar;
figure;s2=pcolor(qxscan,qyscan,EmitScan');s2.FaceColor='interp';s2.EdgeColor='none';xlabel('Qx');ylabel('Qy');title('Emit [pm rad]');colorbar;
figure;s3=pcolor(qxscan,qyscan,BetaX0Scan');s3.FaceColor='interp';s3.EdgeColor='none';xlabel('Qx');ylabel('Qy');title('BetaX0 [m]');colorbar;
figure;s4=pcolor(qxscan,qyscan,BetaY0Scan');s3.FaceColor='interp';s4.EdgeColor='none';xlabel('Qx');ylabel('Qy');title('BetaY0 [m]');colorbar;
    
end
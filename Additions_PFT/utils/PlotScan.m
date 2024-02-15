function PlotScan(rscan,column,label)
%PlotScan Plots results of 2-dim scan of parameters of
%unit cell
Ntheta = rscan.Ntheta;
Nk     = rscan.Nk;

rmat = rscan.data;

figure;
plot(rmat(1:Nk,2),rmat(1:Nk,column));
hold
for i=1:Ntheta-1
    plot(rmat(Nk*i+1:Nk*(i+1),2),rmat(Nk*i+1:Nk*(i+1),column));
end
xlabel('Reverse Bend K [m**-2]');
ylabel(label);

figure;
lines = find(rmat(1:end,2)==rmat(1,2));
plot (rmat(lines,1),rmat(lines,column));
hold;
for j=2:Nk
    krb = rmat(j,2);
    lines = find(rmat(1:end,2)==krb);
    plot (rmat(lines,1),rmat(lines,column));
end
xlabel('Reverse Bend kick [mrad]');
ylabel(label);






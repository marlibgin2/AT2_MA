function [nux, nuy, dpp, maxDnuxP, maxDnuxM, maxDnuyP, maxDnuyM] = TuneEnergyDependence(Rin)

N= 19;
dpp_min = -0.04; dpp_MAX = 0.04; 
dpp = linspace(dpp_min, dpp_MAX, N);


for i=1:length(dpp)
%     tw = gettwiss(Rin, dpp(i));
%     nux(i) = tw.phix(end);
%     nuy(i) = tw.phiy(end);
    [~,lindata]=atlinopt6(Rin,1:length(Rin)+1,'dp',dpp(i));
    nu   = lindata(end).mu(1:2)/2/pi;
    nux(i) = nu(1) * 20;
    nuy(i) = nu(2) * 20;
end

figure(1111); clf; hold on; grid on
plot(dpp, nux-nux(10),'o--');
plot(dpp, nuy-nuy(10),'o--');
figure(1112); clf; hold on; grid on
plot(nux(1:9), nuy(1:9),'b--o','LineWidth',3)
plot(nux(11:19), nuy(11:19),'r--o','LineWidth',3)
axis([48 48.5 14 14.5])

figure(36); hold on;
plot( nuy(1:9)-floor(nuy(10)),'co','LineWidth',3)
plot(nux(10)-floor(nux(10)), nuy(10)-floor(nuy(10)),'y*','LineWidth',3)
plot(nux(11:19)-floor(nux(10)), nuy(11:19)-floor(nuy(10)),'ro','LineWidth',3)

maxDnuxP = max(abs(nux(1:9)  - nux(10)));
maxDnuxM = max(abs(nux(11:19)- nux(10)));
maxDnuyP = max(abs(nuy(1:9)  - nuy(10)));
maxDnuyM = max(abs(nuy(11:19)- nuy(10)));

end
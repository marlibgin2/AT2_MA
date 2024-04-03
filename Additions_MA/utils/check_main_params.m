function RP = check_main_params(ring, per)
%clear THERING
%global THERING
%THERING = ring;

[~,lindata]=atlinopt6(ring,1:length(ring)+1);

if nargin<2
    per=1;
end
%
% RP = ringpars
%
mu   = cat(1,lindata.mu)/2/pi; mux = mu(:,1); muy= mu(:,2);

if 1 == 0
s2   = findcells(ring,'FamName','S2');
s5   = findcells(ring,'FamName','S5');
s4   = findcells(ring,'FamName','S4');
ns4_2   = numel(s4)/1; ns4_1 = numel(s4)/2;
ns2_2   = numel(s2)/1; ns2_1 = numel(s2)/2;
ns5_2   = numel(s5)/1; ns5_1 = numel(s5)/2;

dmux_24 =  (mux(s4(ns4_1))-mux(s2(ns2_1)));
dmuy_24 =  (muy(s4(ns4_1))-muy(s2(ns2_1)));
dmux_45 =  (mux(s5(ns5_1))-mux(s4(ns4_1)));
dmuy_45 =  (muy(s5(ns5_1))-muy(s4(ns4_1)));
dmux_55 =  (mux(s5(ns5_2))-mux(s5(ns5_1)));
dmuy_55 =  (muy(s5(ns5_2))-muy(s5(ns5_1)));
end

% dmux_24 =  (mux(s4(2))-mux(s2(2)))/2/pi;
% dmuy_24 =  (muy(s4(2))-muy(s2(2)))/2/pi;
% dmux_45 =  (mux(s5(2))-mux(s4(2)))/2/pi;
% dmuy_45 =  (muy(s5(2))-muy(s4(2)))/2/pi;
% dmux_55 =  (mux(s5(4))-mux(s5(2)))/2/pi;
% dmuy_55 =  (muy(s5(4))-muy(s5(2)))/2/pi;

AAA  = atsummary_ring(ring);
RP.damping = AAA.damping;
RP.naturalEnergySpread = AAA.naturalEnergySpread;
RP.emix = AAA.naturalEmittance;
RP.xi   = AAA.chromaticity;
RP.nu   = lindata(end).mu(1:2)/2/pi; % AAA.tunes;
RP.U0   = AAA.radiation*1e6; % energy loss in keV
beta = cat(1,lindata.beta);
eta  = cat(2,lindata.Dispersion);
etax = eta(1,:);
bx   = beta(:,1);
by   = beta(:,2);
RP.beta0 = [bx(1) by(1)];
RP.etax  = etax(1);

RP.C    = AAA.circumference;
RP.alfa_C = AAA.compactionFactor;

BBB = findcells(ring,'Class','Bend');
AB  = findcells(ring,'Class','Multipole');
Acell = 0;
for i=1:length(BBB)
    Acell = Acell + ring{BBB(i)}.BendingAngle;
end
for i=1:length(AB)
    Acell = Acell + ring{AB(i)}.PolynomB(1);
end


% tw   = twissring(THERING,0.001,1:length(THERING)+1);
% beta = cat(1,tw.beta);
% bx   = beta(:,1);
% by   = beta(:,2);

RP.period = per; 
RP.TotAng = Acell*RP.period*180/pi;
disp(['emix   = ' num2str(RP.emix*1e12) ' (pm)'])
disp(['  U0   = ' num2str(RP.U0*RP.period) ' (keV)'])
disp(['TotAng = ' num2str(RP.TotAng,9) ' (deg)'])
disp(['TotLen = ' num2str(RP.C*RP.period,9) ' (m)'])

disp(['BETAx,y   = ' num2str(bx(1)) '/'  num2str(by(1))  ' (m)'])

RP.nup = RP.nu*RP.period; RP.xip = RP.xi*RP.period;
disp(['  nu = ' num2str(RP.nup,8) ''])
disp(['  xi = ' num2str(RP.xip,8) ''])
disp([' a_C = ' num2str(RP.alfa_C,8) ''])

if 1 == 0
disp([' dmu_24 = ' num2str(dmux_24,5) ' ' num2str(dmuy_24,5) ' th.: [' num2str(3/7) ' ' num2str(1/7) ' ]'])
disp([' dmu_45 = ' num2str(dmux_45,5) ' ' num2str(dmuy_45,5) ' th.: [' num2str(3/7) ' ' num2str(1/7) ' ]'])
disp([' dmu_55 = ' num2str(dmux_55,5) ' ' num2str(dmuy_55,5) ' th.: [' num2str(3/7) ' ' num2str(1/7) ' ]'])
end

stiff = sqrt(sum((RP.xi./RP.nu).^2));

disp(['  stiffness = ' num2str(stiff) ''])
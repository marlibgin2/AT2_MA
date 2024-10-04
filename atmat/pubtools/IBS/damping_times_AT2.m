function [tx,ty,ts]=damping_times_AT2(ring)

%global THERING GLOBVAL

Ca = 2113.1;        % m^2/GeV^3/s
E=3;   % GeV

L0_tot=0;
for loop=1:length(ring)
   L0_tot=L0_tot+ring{loop}.Length;
end

BENDINDEX = findcells(ring,'BendingAngle');
bend_angle = getcellstruct(ring,'BendingAngle',BENDINDEX);
bend_l = getcellstruct(ring,'Length',BENDINDEX);
rho = bend_l./bend_angle;
max(rho)
sum(rho.*bend_l)./sum(bend_l)
min(rho)
k = getcellstruct(ring,'PolynomB',BENDINDEX,2);

ring = atdisable_6d(ring);
% radiationoff
% cavityoff
dispersion = (findorbit4(ring,0.001,1:length(ring))-findorbit4(ring,-0.001,1:length(ring)))/0.002;

dispersion = dispersion(1,BENDINDEX)';

I2 = sum(bend_l./(rho.^2))
I4x = sum(dispersion./(rho.^3).*(1+2.*(rho.^2).*k).*bend_l)
I4y = 0;

tx = (Ca./L0_tot.*(E.^3).*I2.*(1-I4x./I2)).^-1;
ty = (Ca./L0_tot.*(E.^3).*I2.*(1-I4y./I2)).^-1;
ts = (Ca./L0_tot.*(E.^3).*I2.*(2+(I4x+I4y)./I2)).^-1;
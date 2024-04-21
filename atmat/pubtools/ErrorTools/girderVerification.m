NEWRING = applyGirderError(RING,gI,0,50e-6,50e-6,0.1e-4,0.1e-4,0.1e-3);

s = findspos(NEWRING,1:numel(NEWRING)+1);


RES1 = getcellstruct(NEWRING,'T1',gI);
RES2 = getcellstruct(NEWRING,'T2',gI);
ROT1 = getcellstruct(NEWRING,'R1',gI);
ROT2 = getcellstruct(NEWRING,'R2',gI);
RES1 = cat(2,RES1{:}); RES2 = cat(2,RES2{:});
orb0 = findorbit6(RING,1:numel(RING)+1);
orb = findorbit6(NEWRING,1:numel(NEWRING)+1);
close all

figure; plot(s,orb(1,:)-orb0(1,:),'-b'); hold on; plot(s,orb(3,:)-orb0(3,:),'-r'); xlabel('s [m]'); ylabel('\Deltax,y [m]');  xlabel('s [m]'); title('Closed orbit change')

% Plot translations
figure; 
subplot(2,3,1); plot(s(gI), 1e6*RES1(1,:),'-xb'); hold on; plot(s(gI+1), 1e6*RES2(1,:),'-or'); ylabel('T1(1), T2(1)');  xlabel('s [m]'); title('x [µm]');
subplot(2,3,2); plot(s(gI), 1e6*RES1(3,:),'-xb'); hold on; plot(s(gI+1), 1e6*RES2(3,:),'-or'); ylabel('T1(3), T2(3)');  xlabel('s [m]'); title('y [µm]');
subplot(2,3,3); plot(s(gI), 1e6*RES1(6,:),'-xb'); hold on; plot(s(gI+1), 1e6*RES2(6,:),'-or'); ylabel('T1(6), T2(6)'); xlabel('s [m]'); title('cT [µm]');
subplot(2,3,4); plot(s(gI), RES1(2,:),'-xb'); hold on; plot(s(gI+1), RES2(2,:),'-or'); ylabel('T1(2), T2(2)'); xlabel('s [m]'); title('p_x');
subplot(2,3,5); plot(s(gI), RES1(4,:),'-xb'); hold on; plot(s(gI+1), RES2(4,:),'-or'); ylabel('T1(4), T2(4)'); xlabel('s [m]'); title('p_y');
subplot(2,3,6); plot(s(gI), RES1(5,:),'-xb'); hold on; plot(s(gI+1), RES2(5,:),'-or'); ylabel('T1(5), T2(5)'); xlabel('s [m]'); title('dp/p');

% Plot rotations
phi1 = cellfun(@(x) asin(x(1,3)),ROT1);
phi2 = cellfun(@(x) asin(x(1,3)),ROT2);
figure; phi_ax = plotyy(s(gI), 1e6*phi1, s(gI+1), 1e6*phi2); xlabel('s [m]'); ylabel(phi_ax(1),'R1 \phi [µrad]'); ylabel(phi_ax(2),'R2 \phi [µrad]'); title('Roll angle \phi');

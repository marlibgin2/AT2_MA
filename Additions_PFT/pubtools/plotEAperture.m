function  EApertures=plotEAperture(RING)
% Plots the Elliptical Aperture fields 
% 
nelem = length(RING);
Spos=findspos(RING,1:nelem);
EApertures=zeros(nelem,2);
for i=1:nelem
    if (isfield(RING{i},'EApertures'))
        EApertures(i,1:2) = RING{i}.EApertures;
    else
        EApertures(i,1:2) = nan(1,2);
    end
end
figure;plot(Spos,EApertures(:,1)*1000,'-ob');hold on;xlabel('S[m]');ylabel('Aperture[mm]');
ylim([0,20]);
plot(Spos,EApertures(:,2)*1000,'-sr');
legend('X','Y');grid;
end
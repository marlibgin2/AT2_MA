function  plotRAperture(RING)
% Plots the Rectangular aperture fields 
% 
nelem = length(RING);
Spos=findspos(RING,1:nelem);
RApertures=zeros(nelem,4);
for i=1:nelem
    if (isfield(RING{i},'RApertures'))
        RApertures(i,1:4) = RING{i}.RApertures;
    else
        RApertures(i,1:4) = nan(1,4);
    end
end
figure;plot(Spos,RApertures(:,1)*1000,'-b');hold on;xlabel('S[m]');ylabel('Aperture[mm]');
plot(Spos,RApertures(:,3)*1000,'-r');
plot(Spos,RApertures(:,2)*1000,'-b');
plot(Spos,RApertures(:,4)*1000,'-r');
legend('X','Y');
end
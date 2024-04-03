function rdt = obj_rdt(X, rrr)
%global rrr

qfii    = findcells(rrr,'FamName','QFI');
qfoi    = findcells(rrr,'FamName','QFO');
qfmi   = findcells(rrr,'FamName','QFM');
qfendi = findcells(rrr,'FamName','QFEND');
qdendi = findcells(rrr,'FamName','QDEND');

sdi     = findcells(rrr,'FamName','SD');
sdendi  = findcells(rrr,'FamName','SDEND');
sfmi    = findcells(rrr,'FamName','SFM');
sfoi    = findcells(rrr,'FamName','SFO');
sfii    = findcells(rrr,'FamName','SFI');

dipi   = findcells(rrr,'FamName','DIP');
dipmi  = findcells(rrr,'FamName','DIPm');
for i=1:length(dipi)
    wdip(i) = rrr{dipi(i)}.PolynomB(2)/rrr{dipi(7)}.PolynomB(2);
end
for i=1:length(dipmi)
    wdipm(i) = rrr{dipmi(i)}.PolynomB(2)/rrr{dipmi(7)}.PolynomB(2);
end

VARi = {qfii; qfoi; qfmi; qfendi; qdendi};

for j = 1:5
for i=1:length(VARi{j})
    rrr{VARi{j}(i)}.PolynomB(2) = X(j);
    rrr{VARi{j}(i)}.K           = X(j);
end
end

clear VARi;
VARi = {dipi; dipmi};
W    = {wdip; wdipm};
for j = 1:2
    for i=1:length(VARi{j})
        rrr{VARi{j}(i)}.PolynomB(2) = X(j+5) * W{j}(i);
        rrr{VARi{j}(i)}.K           = X(j+5) * W{j}(i);
    end
end

VARi = {sdi; sdendi; sfmi; sfoi; sfii};
for j = 1:5
for i=1:length(VARi{j})
    rrr{VARi{j}(i)}.PolynomB(3) = X(j+7);
    rrr{VARi{j}(i)}.K           = X(j+7);
end
end

%
% force chromaticity to be 1/1
%
%try
%    rrr = atfitchrom(rrr,[1 1],'SFM','SDEND');
%catch
%    disp('cannot fit chromaticity ... mmachine unstable')
%end

    try
        %RDT = computeRDT_MS(rrr,1:254,'geometric1','geometric2','chromatic','tuneshifts','reuselinear');
        RDT = computeRDT_MS(rrr,1:254,'geometric1','geometric2','reuselinear');
        nfield = length(fieldnames(RDT));

        rdt = 0; fn = fieldnames(RDT);
        for i = 1:nfield
            ff = fn{i};
    
            if  ~strcmpi(ff,'h11001') && ~strcmpi(ff,'h00111') % not the chromaticity! 
                rdt = rdt + abs(RDT(end).(ff));
            end

        end
    catch
        rdt = 1e19;
    end
        
    if isnan(rdt)
        rdt = 1e19;
    end

    disp(['Sum(rdt) = ' num2str(rdt)])

end
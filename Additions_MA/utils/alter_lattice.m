function rrr = alter_lattice(X, rrr)

if ~isempty(X)
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

    oxxoi = findcells(rrr,'FamName','oxxo');
    oxyoi = findcells(rrr,'FamName','oxyo');
    oyyoi = findcells(rrr,'FamName','oyyo');
    
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

    VARi = {oxxoi; oxyoi; oyyoi};

    % zero all the octupoles )in force from 24-11-2023)
    noOCT=1;
    if noOCT==1
        for j=1:3
            for i=1:length(VARi{j})
                rrr{VARi{j}(i)}.PolynomB(4) = 0;
            end
        end
    end
end

end
function MA_RDT = MAcalcRDT(ring)

%
%
%

% [ringdata,lindata]=atlinopt6(ring,1:length(ring)+1);
% beta = cat(1,lindata.beta);
% disp = cat(2,lindata.Dispersion); etax = disp(1,:)';
% bx   = beta(:,1);
% by   = beta(:,2);
S    = findspos(ring,1:length(ring));

[~,AVEBETA,AVEMU,AVEDISP,~,~]=atavedata(ring,0,1:length(ring));
bx=AVEBETA(:,1);
by=AVEBETA(:,2);
etax=AVEDISP(:,1);

for i = 1:length(ring)
    b2L = 0;
    b3L = 0;
    if strcmpi(ring{i}.Class,'Quadrupole') || strcmpi(ring{i}.Class,'Bend')
        b2  =   ring{i}.PolynomB(2) * ring{i}.Length;
        b2L =   b2 * ring{i}.Length;
    elseif strcmpi(ring{i}.Class,'Sextupole')
        b3  =   ring{i}.PolynomB(3) * ring{i}.Length / 2;
        b3L =   b3 * ring{i}.Length; 
    end
    h11001(i) =  1/4 * ((b2L - 2*b3L * etax(i)) * bx(i));
    h00111(i) = -1/4 * ((b2L - 2*b3L * etax(i)) * by(i));

end

9
end
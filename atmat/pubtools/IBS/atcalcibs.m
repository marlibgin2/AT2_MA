function [ex_ibs, ey_ibs, ex0, ey0] = atcalcibs(achro, Ibeam, HCfactor)
%
%
%
verbo    = false;
c        = 2.99793e8;		% m/s
qe       = 1.60218e-19;
h        = 176;
T        = 528/c;
qb       = T * Ibeam / h; % bunch charge (C)
if nargin < 2
    Ibeam = 0.5; % 500 mA
    HCfactor = 5.;
end
ring = achromat2ring(achro);
ring = atdisable_6d(ring);
twissdata  = twissring(ring,0,1:length(ring));
dispersion = (findorbit4(ring,0.001,1:length(ring))-findorbit4(ring,-0.001,1:length(ring)))/0.002;
spos       = findspos(ring,1:length(ring));

ring       = atenable_6d(ring);
BINDEX     = findcells(ring,'PassMethod','BndMPoleSymplectic4RadPass');
QSOINDEX = findcells(ring,'PassMethod','StrMPoleSymplectic4Pass'); % quad+sext+oct ...
for i = 1:numel(QSOINDEX)
    ring{QSOINDEX(i)}.PassMethod =  'StrMPoleSymplectic4RadPass';
end
RADELEMINDEX = sort([BINDEX QSOINDEX]); % radiating elements
[ENV, DP, DL] = ohmienvelopeMA(ring,RADELEMINDEX, 1:length(ring)); % check pb in ohmienvelop with Wig
sigmas = cat(2,ENV.Sigma);
ring   = atdisable_6d(ring);
INJXEMITTANCE = (sigmas(1,1)^2-(dispersion(1,1)*DP)^2)/twissdata(1).beta(1);          % pi m rad
KAPPA         = 0.0125; % betatron  coupling
KAPPAE        = 0.0125; % emittance coupling
YDISPrms      = norm(dispersion(3,:))/sqrt(length(dispersion));  % rms vertical dispersion
XDISPrms      = norm(dispersion(1,:))/sqrt(length(dispersion));  % rms horizontal dispersion
INJSPREAD     = DP;                % Starting energy spread
energye = 3;		% GeV
gam0 = energye/(0.511e-3);
bet0 = sqrt( 1 - gam0.^(-2));

kappa   = KAPPA;		            % betatron-coupling
disprat = (YDISPrms)/(XDISPrms);	% ratio of the dispersion functions

%L = max(spos);			% lenth of the ring
Nemin  = 0.0;
Nemax  = qb / 1.601e-19; % max number of electrons per bunch

betxe = twissdata(1).beta(1);
betye = twissdata(1).beta(2);

[tx,ty,ts]=damping_times_AT2(ring);

taux = (tx)/2;
tauy = (ty)/2;
taul = (ts)/2;

epsnatx = INJXEMITTANCE;
sigenat = INJSPREAD;

siglenat = HCfactor * DL; % bunch lengthening from HC

epslnat = siglenat * sigenat;
betle = siglenat / sigenat;

Nel=[];epsxel=[];epsyel=[];deltapel=[];siglel=[];

epsxe = INJXEMITTANCE/(1+KAPPAE);
epsye = INJXEMITTANCE*KAPPAE/(1+KAPPAE);
deltape = INJSPREAD;
epsle = (deltape.^2.0) * betle;

time = 0.0;

disp('Initial Values');

sigxe = sqrt( epsxe*betxe);
sigye = sqrt( epsye*betye);
sigle = sqrt( epsle*betle);
deltape = sqrt( epsle/betle);
ex0 = epsxe; ey0 = epsye;
fprintf('Eps_x = %.3g pm.rad, Eps_y = %.3g pm.rad, Delta p/p = %.3e\n', epsxe*1e12, epsye*1e12, deltape);
fprintf('Sig_x = %.3g um, Sig_y = %.3g um, Sig_z = %.3g mm\n\n', sigxe*1e6, sigye*1e6, sigle*1000);

Ne = Nemax;
while 1


    if (Ne > 0.0)
        [ibsle,ibsxe,ibsye]=ibs_bane(ring,Ne,sigle,deltape,epsxe,epsye);
    else
        ibsle=1e100;
        ibsxe=1e100;
        ibsye=1e100;
    end


    % Emittanzaenderung ausrechnen

    depsxe = -1/taux * ( epsxe - 1/(1+kappa)*epsnatx) + 1/ibsxe * epsxe;
    depsye = -1/tauy * ( epsye - (2*kappa/(1+kappa)*epsxe+disprat.^2*epsnatx)) + 1/ibsye * epsye + 1/ibsxe * (disprat.^2) * epsxe;
    depsle = -1/taul * ( epsle - epslnat) + 1/ibsle * epsle;

    timestep = 1/( 1/taul + 1/ibsle)/4;

    % Aendere Emittanz in linearer Naeherung

    epsxe = epsxe + depsxe * timestep;
    epsye = epsye + depsye * timestep;
    epsle = epsle + depsle * timestep;

    time = time + timestep;

    sigxe = sqrt( epsxe*betxe);
    sigye = sqrt( epsye*betye);
    sigle = sqrt( epsle*betle);
    deltape = sqrt( epsle/betle);
    if verbo
        fprintf( 't = %.3g s, Eps_x = %.4g nm rad, Eps_y = %.4g nm rad, Delta p/p = %.3e\n', time, epsxe*1e9, epsye*1e9, deltape);
    end
    if ( ( abs( (depsxe*timestep)/epsxe) < 1.0e-4) & ...
            ( abs( (depsye*timestep)/epsye) < 1.0e-4) & ...
            ( abs( (depsle*timestep)/epsle) < 1.0e-4))
        break;
    end

    if ( time > 1.0)
        break;	% maximal 1 s
    end
end
disp ('Final values')
fprintf('Ne = %.3g \n',Ne)
fprintf('Eps_x = %.3g pm.rad, Eps_y = %.3g pm.rad, Delta p/p = %.3e\n', epsxe*1e12, epsye*1e12, deltape);
fprintf('Sig_x = %.3g um, Sig_y = %.3g um, Sig_z = %.3g mm\n\n', sigxe*1e6, sigye*1e6, sigle*1000);
ex_ibs = epsxe; ey_ibs=epsye;
end
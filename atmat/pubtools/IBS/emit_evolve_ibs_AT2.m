function emit_evolve_ibs_AT2(ring)
% Iterative Berechnung der Emittanzentwicklung bei konstanter Energie
%
%	Beruecksichtigt werden Strahlungsdaempfung, transversale Kopplung,
%	Intrabeam-Scattering (wird mit ibs_bane.m berechnet)
%
%	(w) 1998 Michael Gentner, Phys.Inst.Univ.Bonn
%           2000-2010 Christoph Steier, Lawrence Berkeley National Laboratory

c = 2.99793e8;		% m/s
qe = 1.60218e-19;   % C
me = 9.1081e-31;    % kg
IA = 17.0e3;		% A
re = 2.818e-15;		% m

%global THERING GLOBVAL

%%% maxivbare_SextOct;
% maxiv_5IDs;

ring = achromat2ring(ring); 
ring = atdisable_6d(ring);
twissdata  = twissring(ring,0,1:length(ring));
dispersion = (findorbit4(ring,0.001,1:length(ring))-findorbit4(ring,-0.001,1:length(ring)))/0.002;
spos       = findspos(ring,1:length(ring));

ring       = atenable_6d(ring);
BENDINDEX     = findcells(ring,'PassMethod','BndMPoleSymplectic4RadPass');
QUADSEXTINDEX = findcells(ring,'PassMethod','StrMPoleSymplectic4Pass');
for i = 1:numel(QUADSEXTINDEX)
   ring{QUADSEXTINDEX(i)}.PassMethod =  'StrMPoleSymplectic4RadPass';
end
RADELEMINDEX = sort([BENDINDEX QUADSEXTINDEX]);
[ENV, DP, DL] = ohmienvelopeMA(ring,RADELEMINDEX, 1:length(ring)); % check pb in ohmienvelop with Wig
sigmas = cat(2,ENV.Sigma);

ring = atdisable_6d(ring);
INJXEMITTANCE=(sigmas(1,1)^2-(dispersion(1,1)*DP)^2)/twissdata(1).beta(1);          % pi m rad
KAPPA=0.0125; %Bare lattice value0.0125;          % Betatron-Kopplung
% for maxiv_5IDs lattice KAPPA=0.014;
KAPPAE=0.0125;%Bare lattice value0.0125;       	  % Emittanz-Kopplung
YDISPrms=norm(dispersion(3,:))/sqrt(length(dispersion));;                    % rms vertical dispersion
XDISPrms=norm(dispersion(1,:))/sqrt(length(dispersion));                % rms horizontal dispersion
INJSPREAD=DP;                % Starting energy spread

energye = 3;		% GeV
gam0 = energye/(0.511e-3);
bet0 = sqrt( 1 - gam0.^(-2));

kappa = KAPPA;		% Betatron-Kopplung
disprat = (YDISPrms)/(XDISPrms);	% Verhaeltnis der Dispersionsfunktionen

L = max(spos);			% Umfang Max-5

Nenom=1.5615e10; % Equivelant to 250 mA
Nemin=0.0;
Nemax=2*Nenom; % Here it is 500 mA
Nestep=Nenom/5;

betxe = twissdata(1).beta(1);
betye = twissdata(1).beta(2);
% radiationoff
% cavityoff
% twissdata=twissring(THERING,0,1:length(THERING));
% dispersion = (findorbit4(THERING,0.001,1:length(THERING))-findorbit4(THERING,-0.001,1:length(THERING)))/0.002;
% spos=findspos(THERING,1:length(THERING));

% radiationon
% cavityon
% BENDINDEX = findcells(THERING,'PassMethod','BndMPoleSymplectic4RadPass');
% QUADSEXTINDEX = findcells(THERING,'PassMethod','StrMPoleSymplectic4RadPass');
% RADELEMINDEX = sort([BENDINDEX QUADSEXTINDEX]);
% [ENV, DP, DL] = ohmienvelope(THERING,RADELEMINDEX, 1:length(THERING));
% sigmas = cat(2,ENV.Sigma);
% radiationoff;
% cavityoff;


[tx,ty,ts]=damping_times_AT2(ring);

taux = (tx)/2;
tauy = (ty)/2;
taul = (ts)/2;

epsnatx = INJXEMITTANCE;
sigenat = INJSPREAD;
siglenat = 5.545*DL; % Hamed: LC Lengthening with bare lattice to reach 56mm is 5.545
siglenat = 1*DL; % Hamed: LC Lengthening with bare lattice to reach 56mm is 5.545
%siglenat = DL; % MA bare lattice no HC lengthening
epslnat = siglenat * sigenat;
betle = siglenat / sigenat;

loop=1;
figure;

Nel=[];epsxel=[];epsyel=[];deltapel=[];siglel=[];

for Ne = Nemin:Nestep:Nemax

    % Startwerte

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

    fprintf('Eps_x = %.3g nm rad, Eps_y = %.3g nm rad, Delta p/p = %.3e\n', epsxe*1e9, epsye*1e9, deltape);
    fprintf('Sig_x = %.3g mm, Sig_y = %.3g mm, Sig_z = %.3g mm\n\n', sigxe*1000, sigye*1000, sigle*1000);

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

        fprintf( 't = %.3g s, Eps_x = %.4g nm rad, Eps_y = %.4g nm rad, Delta p/p = %.3e\n', time, epsxe*1e9, epsye*1e9, deltape);

        if ( ( abs( (depsxe*timestep)/epsxe) < 1.0e-4) & ...
                ( abs( (depsye*timestep)/epsye) < 1.0e-4) & ...
                ( abs( (depsle*timestep)/epsle) < 1.0e-4))
            break;
        end

        if ( time > 1.0)
            break;	% maximal 1 s
        end
    end



    fprintf( '%e\t%e\t%e\t%e\t%e\n', Ne, epsxe, epsye, deltape, sigle);

    Nel(loop)=Ne;epsxel(loop)=epsxe;epsyel(loop)=epsye;deltapel(loop)=deltape;siglel(loop)=sigle;
    loop=loop+1;

    if loop>2
        subplot(2,2,1);
        plot(Nel*(176*c *qe/528),epsxel,'.-','LineWidth',2);
        set(gca,'FontSize',12);
        xlabel('Beam Current [A]','FontSize',14);
        ylabel('\epsilon_x [m.rad]','FontSize',14);
        grid on
        subplot(2,2,2);
        plot(Nel*(176*c *qe/528),epsyel,'.-','LineWidth',2);
        set(gca,'FontSize',12);
        xlabel('Beam Current [A]','FontSize',14);
        ylabel('\epsilon_y [m.rad]','FontSize',14);
        grid on
        subplot(2,2,3);
        plot(Nel*(176*c *qe/528),deltapel,'.-','LineWidth',2);
        set(gca,'FontSize',12);
        xlabel('Beam Current [A]','FontSize',14);
        ylabel('\Delta p/p','FontSize',14);
        grid on
        drawnow;
        subplot(2,2,4);
        plot(Nel*(176*c *qe/528),siglel,'.-','LineWidth',2);
        set(gca,'FontSize',12);
        xlabel('Beam Current [A]','FontSize',14);
        ylabel('\sigma_{z} [m]','FontSize',14);
        grid on
        drawnow;
        pause(0.1);
    end
end


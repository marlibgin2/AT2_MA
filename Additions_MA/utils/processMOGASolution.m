clear EM hcmi vcmi sVCM sHCM HCM
RI     = m4U_240314_b01_02_04_03;  
load MOGA_24032024_192158.mat; 
RI     = alter_m4U_lattice(X(15,:), RI, 'B'); clear RI; 
Rnu    = atfittune(Rout,0,[55.40/20, 17.10/20],'Q1','Q2');
Rnuch  = atfitchrom(Rnu,[1/20, 1/20],'S1','S2');
clear RI Rnu
A  = Rnuch; clear Rnuch
A  = atenable_6d(A);
TL = achromat2ring(A);
TL = atenable_6d(TL); 
CAVi = findcells(TL,'FamName','CAV');
for ii=1:6; TL{CAVi(ii)}.Voltage=1.8e6/6; end
hcmi = findcells(TL,'iscorH'); 
vcmi = findcells(TL,'iscorV');
sHCM = findspos(TL,hcmi);
sVCM = findspos(TL,vcmi);

r0     = 3e-3;
nsteps = 20; 
nturns = 2000;
dp     = 0.0;
res    = 0.25e-3;
alpha  = 1.1;
     
for js = 1:10
    disp(['running seed' num2str(js)])
    EM{js}   = errormodel_Gradients_reducedMisalignment; %errormodel_onlyGradients; %    
    okappa = 0;
    while okappa == 0
    try
        TLe{js}  = applyErrorModel(TL,EM{js});
        okappa = 1;
    catch
        okappa = 0;
    end
    end

     TLec{js} = atcorrectorbit(TLe{js},[],[],[],[], ...
         [140 120; 160 140; 180 160; ones(10,1)*[200 180]],...
         [true true],0.75,[],[],[0.38, 0.38]*1e-3);

    % TLec{js} = TLe{js};
    [o, s] = getorbit_everywhere_maxiv_at2(TLec{js});
    figure(1); clf; hold on
    plot(s, 1e6*o(1,:),'b');
    plot(s, 1e6*o(3,:),'r');
    axis([0 528 -300 300]); xlabel('S (m)'); ylabel('CM (urad)')
    hcmi = findcells(TLec{js},'iscorH'); vcmi = findcells(TLec{js},'iscorV');
    for i=1:length(hcmi); hcm(i) = TLec{js}{hcmi(i)}.KickAngle(1); end
    for i=1:length(vcmi); vcm(i) = TLec{js}{vcmi(i)}.KickAngle(2); end
    figure(2); clf; hold on
    stem(sHCM, hcm*1e6,'b','linewidth',3,'marker','none');
    stem(sVCM, vcm*1e6,'r','linewidth',3,'marker','none');
    axis([0 528 -350 350]); xlabel('S (m)'); ylabel('CM (urad)')
    O{js} = o; HCM{js}=hcm; VMC{js}=vcm;

    DA{js} = modelDA_sim_par(TLec{js}, r0, nsteps, nturns, dp, res, alpha);
end

figure(335); clf; hold on
for js = 1:10
    plot(1e3*DA{js}(:,1),1e3*DA{js}(:,2))
end
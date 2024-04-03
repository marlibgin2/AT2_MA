if 1 == 1
Rout   = ExploreMOGASolution_m4U('MOGA_24032024_192158.mat',15);
Rnu    = atfittune(Rout,0,[55.40/20, 17.10/20],'Q1','Q2');
Rnuch  = atfitchrom(Rnu,[1/20, 1/20],'S1','S2');

TL = achromat2ring(Rnuch); 
% tic; [x0pos, y0pos, nuxpos, nuypos, diffuvec, WA] = DA_NU_FMA_par(TL,'fma',2000,200/4,128/4,12e-3,7e-3,'newcalc','cluster cores',128,'graphic_plots'); toc
% figure(36); axis([0 0.5 0 0.5])

EM = errormodel_onlyGradients;
TLe = applyErrorModel(TL,EM);
pause()
tic; [x0pos, y0pos, nuxpos, nuypos, diffuvec, WA] = DA_NU_FMA_par(TLe,'fma',2000,200/4,128/4,12e-3,7e-3,'newcalc','cluster cores',128,'graphic_plots'); toc
figure(36); axis([0 0.5 0 0.5])

end
r0     = 3e-3;
       nsteps = 20; %61;
        nturns = 2000;
        dp     = 0.0;
        res    = 0.25e-3;
        alpha  = 1.1;
        
        %seed_case = {EM1, EM2, EM3, EM4}; 
        %TLe_case  = {TLe_s1, TLe_s2, TLe_s3, TLe_s4}; 

        if 1 == 1
        %clear DA;  EMseed={}; TLeseed = {};
        for i = 1:10
            EMseed{i}  = errormodel_onlyGradients; 
            pause(2)
            TLeseed{i} = applyErrorModel(TL, EMseed{i});    
            disp(['examining seed n. ' num2str(i)])
            pause(2)
            DA{i} = modelDA_sim_par(TLeseed{i}, r0, nsteps, nturns, dp, res, alpha);
        end
        end

        %clear DA;  EMseed={}; TLeseed = {};
        for i = 11:20
            EMseed{i}  = errormodel_Gradients_reducedMisalignment; 
            pause(2)
            TLeseed{i} = applyErrorModel(TL, EMseed{i}); 
            TLeseed{i} = atenable_6d(TLeseed{i});
            IS_6D = check_6d(TLeseed{i});
            TLeseedc{i} = atcorrectorbit(TLeseed{i},[],[],[],[],[140 120; 160 140; 180 160; ones(10,1)*[200 180]],[true true],0.75,[],[],[0.38, 0.38]*1e-3);
            disp(['examining seed n. ' num2str(i)])
            pause(2)
            DA{i} = modelDA_sim_par(TLeseedc{i}, r0, nsteps, nturns, dp, res, alpha);
        end




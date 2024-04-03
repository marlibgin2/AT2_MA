

%
%
% change octupoles
OXXOval = [0 -500 -900 -1000 -1100 -1500 -1100 -1000  ];
OXYOval = [0  0    0    0     0     0     0     0     ];
OYYOval = [0  0    0    0     0     0    +300  +300   ];

for io = 1:length(OXXOval)
    for j = 1:2
        RING_matched_optconstr__fitC{oxxoi(j)}.PolynomB(4)=OXXOval(io);
        RING_matched_optconstr__fitC{oyyoi(j)}.PolynomB(4)=OYYOval(io); 
        RING_matched_optconstr__fitC{oxyoi(j)}.PolynomB(4)=OXYOval(io);
    end
    Rout = RING_matched_optconstr__fitC; 
    name_out = ['DA_' num2str(OXXOval(io)) '_' num2str(OXYOval(io)) '_' num2str(OYYOval(io)) '.mat'];
    generaXY_DAgrid;
    DA = modelDA_sim_par( Rout, 12e-3, 53, 2000, 0, 0.25e-3, 1.1);
    [Area, ~] = calcDA_Area(DA);
    figure(333); hold on; plot(DA(:,1),DA(:,2),'r--o','linewidth',3)
    text(0.004, 6e-3, ['A = ' num2str(Area*1e6) '(mm^{2})' ],'color','w')
    [nux, nuy, dpp, maxDnuxP, maxDnuxM, maxDnuyP, maxDnuyM] = TuneEnergyDependence(Rout);
    maxDnux = maxDnuxM+maxDnuxP;
    maxDnuy = maxDnuyM+maxDnuyP;
    text(0.004, 5.5e-3, ['max{\delta\nu}_x = ' num2str(maxDnux)],'color','w')
    text(0.004, 5.0e-3, ['max{\delta\nu}_y = ' num2str(maxDnuy)],'color','w')
    save(name_out,'mDA','DA','Area','maxDnux','maxDnuy');
    disp(['case ' name_out '.... pausing 10sec']); pause(10)
end

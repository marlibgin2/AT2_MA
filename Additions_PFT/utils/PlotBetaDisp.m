function PlotBetaDisp(lindata, stitle)
%PlotBetaDisp Plots Beta function and dispersion stored in lattica data
%structure
%   Detailed explanation goes here

beta=cat(1,lindata.beta);
betax=beta(:,1);
betay=beta(:,2);
Dispersion=cat(2,lindata.Dispersion);
Disp=Dispersion(1,:)';
SPos=cat(1,lindata.SPos);
figure;
yyaxis left; ylabel('BetaXY (m)');
plot(SPos,betax,'b-');hold on; plot(SPos,betay,'r-');
yyaxis right;  ylabel('Disp (m)');
plot(SPos,Disp,'g-');
title(stitle);
xlabel('S(m)');
legend('BetaX','BetaY','Disp')

end


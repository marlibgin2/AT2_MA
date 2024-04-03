function [C,Ceq] = CONfcn_m4U(X, rrr)
verbo = true; 

% -----------------------------------------
% alter the lattice according to parameters
% -----------------------------------------

rrr  = alter_m4U_lattice(X, rrr, 'B');
s = findspos(rrr,1:length(rrr)+1); se = s(end);
% guess the periodicity
P = round(528/se); 

% ----------------------------
% force chromaticity to be 1/1
% ----------------------------
changechro=1;
if changechro == 1
try
    rrr = atfitchrom(rrr,[1/P 1/P],'S1','S2'); 
    if verbo; disp('CON-fitting chromaticity to [1,1]'); end
%    error('cannot fit chromaticity ... machine unstable')
catch
    if verbo; disp('cannot fit chromaticity ... machine unstable'); end
end
end


try

[twiss, tune, chrom]  = twissring(rrr,0.0, 1:length(rrr), 'chrom', 1e-7);
eta   = cat(2,twiss.Dispersion); etax=eta(1,:);
beta  = cat(1,twiss.beta); betax = beta(:,1); betay = beta(:,2); 
hi    = abs(etax(1));
bxi   = abs(betax(1));
byi   = abs(betay(1));
bxmax = max(abs(betax));
bymax = max(abs(betay));
mu=cat(1,twiss.mu)/2/pi; mux = mu(:,1); muy = mu(:,2); 
s2   = findcells(rrr,'FamName','S2');
s5   = findcells(rrr,'FamName','S5');
s4   = findcells(rrr,'FamName','S4');
% % % ns4_2   = numel(s4)/1; ns4_1 = numel(s4)/2;
% % % ns2_2   = numel(s2)/1; ns2_1 = numel(s2)/2;
% % % ns5_2   = numel(s5)/1; ns5_1 = numel(s5)/2;
% % % dmux_24 =  (mux(s4(ns4_1))-mux(s2(ns2_1)));
% % % dmuy_24 =  (muy(s4(ns4_1))-muy(s2(ns2_1)));
% % % dmux_45 =  (mux(s5(ns5_1))-mux(s4(ns4_1)));
% % % dmuy_45 =  (muy(s5(ns5_1))-muy(s4(ns4_1)));
% % % dmux_55 =  (mux(s5(ns5_2))-mux(s5(ns5_1)));
% % % dmuy_55 =  (muy(s5(ns5_2))-muy(s5(ns5_1)));

Xix   = chrom(1); Xiy = chrom(2);
Nux   = tune(1);  Nuy = tune(2); 
% % % Rnu   = Nux / Nuy; 
if isnan(Xix*Xiy) || isnan(Nux*Nuy)
     error('unstable machine -- retry')
end
if verbo; disp('CON-try OK'); end

C(1)   = (hi-100e-6)/10e-6; % C<=0 
C(2)   = (bxi-16)/1;%(bxi-10)/1;
C(3)   = (byi-4.25)/1;%(byi-10)/1;
C(4)   = (bxmax-20)/1;
C(5)   = (bymax-20)/1; 
C(6)   = (mod(nux,1)-0.48)/0.01; 
C(7)   = (mod(nuy,1)-0.48)/0.01; 

C(8)   = -1; %(9.2-bxi)/1; %abs(dphix1 - 3/7)/(3/7)-0.4; % dphix wrt 3/7 less than 30% 
C(9)   = -1; %abs(dmux_24 - 3/7)/(3/7)-0.2; % dmux less than 40%
C(10)  = -1; %abs(dmux_45 - 3/7)/(3/7)-0.2; % dmux less than 40% 
C(11)  = -1; %abs(dmux_55 - 3/7)/(3/7)-0.1; % dmux less than 40%
C(12)  = -1; %abs(dmuy_24 - 1/7)/(1/7)-0.4; % dmux less than 40%
C(13)  = -1; %abs(dmuy_45 - 1/7)/(1/7)-0.3; % dmux less than 40% 
C(14)  = -1; %abs(dmuy_55 - 1/7)/(1/7)-0.2; % dmux less than 40%

% C(15)  = -1; %(abs(Rnu) - 3)/3 - 0.6; %0.4 

Ceq(1) = 0; 
Ceq(2) = 0; 
Ceq(3) = 0; 
Ceq(4) = 0; 
Ceq(5) = 0; 
Ceq(6) = 0; 
Ceq(7) = 0; 
Ceq(8) = 0; 
Ceq(9) = 0; 
Ceq(10) = 0; 
Ceq(11) = 0; 
Ceq(12) = 0;
Ceq(13) = 0; 
Ceq(14) = 0; 
% Ceq(6) = abs(xix - 1)/1;
% Ceq(7) = abs(xiy - 1)/1;

catch

if verbo; disp('CON-catch ...'); end

C(1) = 1e19;
C(2) = 1e19;
C(3) = 1e19;
C(4) = 1e19;
C(5) = 1e19;
C(6) = 1e19;
C(7) = 1e19;
C(8) = 1e19;
C(9) = 1e19;
C(10) = 1e19;
C(11) = 1e19;
C(12) = 1e19; 
C(13) = 1e19;
C(14) = 1e19;

Ceq(1) = 0;
Ceq(2) = 0; 
Ceq(3) = 0; 
Ceq(4) = 0; 
Ceq(5) = 0; 
Ceq(6) = 0; 
Ceq(7) = 0; 
Ceq(8) = 0; 
Ceq(9) = 0; 
Ceq(10) = 0; 
Ceq(11) = 0; 
Ceq(12) = 0; 
Ceq(13) = 0; 
Ceq(14) = 0; 

end
end

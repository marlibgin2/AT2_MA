function [C,Ceq] = CONfcn(X, rrr)
verbo = false; 

% -----------------------------------------
% alter the lattice according to parameters
% -----------------------------------------
rrr  = alter_lattice(X, rrr);
[twiss, tune, chrom]  = twissring(rrr,0.0, 1:length(rrr)+1, 'chrom', 1e-7);

% % % if ~isempty(X)
% % %     qfii    = findcells(rrr,'FamName','QFI');
% % %     qfoi    = findcells(rrr,'FamName','QFO');
% % %     qfmi   = findcells(rrr,'FamName','QFM');
% % %     qfendi = findcells(rrr,'FamName','QFEND');
% % %     qdendi = findcells(rrr,'FamName','QDEND');
% % % 
% % %     sdi     = findcells(rrr,'FamName','SD');
% % %     sdendi  = findcells(rrr,'FamName','SDEND');
% % %     sfmi    = findcells(rrr,'FamName','SFM');
% % %     sfoi    = findcells(rrr,'FamName','SFO');
% % %     sfii    = findcells(rrr,'FamName','SFI');
% % % 
% % %     dipi   = findcells(rrr,'FamName','DIP');
% % %     dipmi  = findcells(rrr,'FamName','DIPm');
% % %     for i=1:length(dipi)
% % %         wdip(i) = rrr{dipi(i)}.PolynomB(2)/rrr{dipi(7)}.PolynomB(2);
% % %     end
% % %     for i=1:length(dipmi)
% % %         wdipm(i) = rrr{dipmi(i)}.PolynomB(2)/rrr{dipmi(7)}.PolynomB(2);
% % %     end
% % % 
% % %     VARi = {qfii; qfoi; qfmi; qfendi; qdendi};
% % % 
% % %     for j = 1:5
% % %         for i=1:length(VARi{j})
% % %             rrr{VARi{j}(i)}.PolynomB(2) = X(j);
% % %             rrr{VARi{j}(i)}.K           = X(j);
% % %         end
% % %     end
% % % 
% % %     clear VARi;
% % %     VARi = {dipi; dipmi};
% % %     W    = {wdip; wdipm};
% % %     for j = 1:2
% % %         for i=1:length(VARi{j})
% % %             rrr{VARi{j}(i)}.PolynomB(2) = X(j+5) * W{j}(i);
% % %             rrr{VARi{j}(i)}.K           = X(j+5) * W{j}(i);
% % %         end
% % %     end
% % % 
% % %     VARi = {sdi; sdendi; sfmi; sfoi; sfii};
% % %     for j = 1:5
% % %         for i=1:length(VARi{j})
% % %             rrr{VARi{j}(i)}.PolynomB(3) = X(j+7);
% % %             rrr{VARi{j}(i)}.K           = X(j+7);
% % %         end
% % %     end
% % % 
% % % end

%
% force chromaticity to be 1/1
%
changechro=1;
if changechro == 1
try
    rrr = atfitchrom(rrr,[1 1],'SFM','SDEND'); 
    if verbo; disp('CON fitting chromaticity to [1,1]'); end
%    error('cannot fit chromaticity ... machine unstable')
catch
    if verbo; disp('cannot fit chromaticity ... machine unstable'); end
end
end


try

    
% % % twiss = gettwiss(rrr, 0.0);
% % % hi    = abs(twiss.etax(1));
% % % bxi   = abs(twiss.betax(1));
% % % byi   = abs(twiss.betay(1));
% % % bxmax = max(abs(twiss.betax));
% % % bymax = max(abs(twiss.betay));

[twiss, tune, chrom]  = twissring(rrr,0.0, 1:length(rrr), 'chrom', 1e-7);
eta  = cat(2,twiss.Dispersion); etax=eta(1,:);
beta = cat(1,twiss.beta); betax = beta(:,1); betay = beta(:,2); 
hi    = abs(etax(1));
bxi   = abs(betax(1));
byi   = abs(betay(1));
bxmax = max(abs(betax));
bymax = max(abs(betay));
mu=cat(1,twiss.mu)/2/pi; phix = mu(:,1); phiy = mu(:,2); 
dphix1 =  phix(96)-phix(66); % phase advance between U1 and U2
dphix2 = phix(126)-phix(96); % phase advance between U2 and U3
dphiy1 =  phiy(96)-phiy(66); % phase advance between U1 and U2
dphiy2 = phiy(126)-phiy(96); % phase advance between U2 and U3

xix   = chrom(1); xiy = chrom(2);
nux   = tune(1);  nuy = tune(2); 
Rnu   = nux / nuy; 
if isnan(xix*xiy) || isnan(nux*nuy)
     error('unstable machine -- retry')
end
if verbo; disp('try OK'); end

C(1)   = (hi-100e-6)/10e-6; % C<=0 
C(2)   = (bxi-16)/1;%(bxi-10)/1;
C(3)   = (byi-10)/1;
C(4)   = (bxmax-20)/1;
C(5)   = (bymax-20)/1; 
C(6)   = (mod(nux,1)-0.48)/0.01; 
C(7)   = (mod(nuy,1)-0.48)/0.01; 

C(8)   = abs(dphix1 - 3/7)/(3/7)-0.4; % dphix wrt 3/7 less than 30% 
C(9)   = abs(dphix2 - 3/7)/(3/7)-0.4; % dphix wrt 3/7 less than 30% 
C(10)  = abs(dphiy1 - 1/7)/(1/7)-0.4; % dphix wrt 3/7 less than 30% 
C(11)  = abs(dphiy2 - 1/7)/(1/7)-0.4; % dphix wrt 3/7 less than 30% 
C(12)  = (abs(Rnu) - 3)/3 - 0.6; %0.4 

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
% Ceq(6) = abs(xix - 1)/1;
% Ceq(7) = abs(xiy - 1)/1;

catch

if verbo; disp('catch ...'); end

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

end
end

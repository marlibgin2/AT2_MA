function [C,Ceq] = CONfcn2(X, rrr)

% -----------------------------------------
% alter the lattice according to parameters
% -----------------------------------------
ralt  = alter_lattice(X, rrr); 
[twiss, tune, chrom]  = twissring(ralt,0.0, 1:length(ralt)+1, 'chrom', 1e-7);
% if ~anynan(tune)&&~anynan(chrom)
%     rrr = ralt;
% end

% ----------------------------
% force chromaticity to be 1/1
% ----------------------------
changechro=1;
if changechro == 1
    try
        rtemp = atfitchrom(ralt,[1 1],'SFM','SDEND');
        disp('CON fitting chromaticity to [1,1]')
        [twiss, tune, chrom]  = twissring(rtemp,0.0, 1:length(rtemp)+1, 'chrom', 1e-7);
        ralt = rtemp; clear rtemp;
    catch
        disp('cannot fit chromaticity ... machine unstable')
    end
end


try
% % % twiss = gettwiss(rrr, 0.0);
% % % hi    = abs(twiss.etax(1));
% % % bxi   = abs(twiss.betax(1));
% % % byi   = abs(twiss.betay(1));
% % % bxmax = max(abs(twiss.betax));
% % % bymax = max(abs(twiss.betay));

[twiss, tune, chrom]  = twissring(ralt,0.0, 1:length(ralt)+1, 'chrom', 1e-7);
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
if isnan(xix*xiy) || isnan(nux*nuy)
     error('unstable machine -- retry')
end
disp('CON: tune/chrom are real numbers - try OK')

C(1)   = (hi-100e-6)/10e-6; % C<=0 
C(2)   = (bxi-10)/1;
C(3)   = (byi-10)/1;
C(4)   = (bxmax-20)/1;
C(5)   = (bymax-20)/1; 
C(6)   = (mod(nux,1)-0.49)/0.01; 
C(7)   = (mod(nuy,1)-0.49)/0.01; 
% C(8)   = abs(dphix1 - 3/7)/(3/7)-0.6; % dphix wrt 3/7 less than 30% 
% C(9)   = abs(dphix2 - 3/7)/(3/7)-0.6; % dphix wrt 3/7 less than 30% 
% C(10)  = abs(dphiy1 - 1/7)/(1/7)-0.6; % dphix wrt 3/7 less than 30% 
% C(11)  = abs(dphiy2 - 1/7)/(1/7)-0.6; % dphix wrt 3/7 less than 30% 


Ceq(1) = 0; 
Ceq(2) = 0; 
Ceq(3) = 0; 
Ceq(4) = 0; 
Ceq(5) = 0; 
Ceq(6) = 0; 
Ceq(7) = 0; 
% Ceq(8) = 0; 
% Ceq(9) = 0; 
% Ceq(10) = 0; 
% Ceq(11) = 0; 
% Ceq(6) = abs(xix - 1)/1;
% Ceq(7) = abs(xiy - 1)/1;

catch

disp('CON: tune/chrom are Nan - catch ...')

C(1) = 1e19;
C(2) = 1e19;
C(3) = 1e19;
C(4) = 1e19;
C(5) = 1e19;
C(6) = 1e19;
C(7) = 1e19;
% C(8) = 1e19;
% C(9) = 1e19;
% C(10) = 1e19;
% C(11) = 1e19;

Ceq(1) = 0;
Ceq(2) = 0; 
Ceq(3) = 0; 
Ceq(4) = 0; 
Ceq(5) = 0; 
Ceq(6) = 0; 
Ceq(7) = 0; 
% Ceq(8) = 0; 
% Ceq(9) = 0; 
% Ceq(10) = 0; 
% Ceq(11) = 0; 

end
end


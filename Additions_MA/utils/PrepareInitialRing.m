function Rout = PrepareInitialRing(X0)
%
% some MOGA starting points
%
% MOGA_30112023_022754.mat   Sol. 18
% load MOGA_30112023_022754.mat
% X0=X; clear X 
% use X0 to define your initial lattice >>> Rout
%


global THERING 

% if ~isfile(maxmed_7BA_blank_AT2_simple.mat)
%if ~isfile('maxivu_A_blank_AT2.mat')
if ~isfile('maxivu_B_blank_AT2.mat')
    return
%    medmax_7BA_2_1_1_AT2_simple
    % maxivu_A_blank_AT2
    %%%%% maxivU_B_blank_AT2  %%%% file not found 6/3/2024
    RING1  = THERING; clear THERING
%    save('maxivu_A_blank_AT2.mat','RING1')
    %%%% save('maxivu_B_blank_AT2.mat','RING1') %%%% commented out 6/3/2024
else
%%    load 'maxivu_A_blank_AT2.mat'
    load 'maxivu_B_blank_AT2.mat'
    
end
%Rtmp = RING1; clear RING1
Rtmp = RINGb; clear RINGb

% -------------------------------------------------------------------------
% vector X contains the parameters to be changed to generate a new lattice
% -------------------------------------------------------------------------

qfii   = findcells(Rtmp,'FamName','QFI');
qfoi   = findcells(Rtmp,'FamName','QFO');
qfmi   = findcells(Rtmp,'FamName','QFM');
qfendi = findcells(Rtmp,'FamName','QFEND');
qdendi = findcells(Rtmp,'FamName','QDEND');

sdi     = findcells(Rtmp,'FamName','SD');
sdendi  = findcells(Rtmp,'FamName','SDEND');
sfmi    = findcells(Rtmp,'FamName','SFM');
sfoi    = findcells(Rtmp,'FamName','SFO');
sfii    = findcells(Rtmp,'FamName','SFI');

dipi   = findcells(Rtmp,'FamName','DIP');
dipmi  = findcells(Rtmp,'FamName','DIPm');

oxxoi = findcells(Rtmp,'FamName','oxxo');
oxyoi = findcells(Rtmp,'FamName','oxyo');
oyyoi = findcells(Rtmp,'FamName','oyyo');

% zero all the octupoles )in force from 24-11-2023)
VARi = {oxxoi; oxyoi; oyyoi};
noOCT=1;
if noOCT==1
    for j=1:3
        for i=1:length(VARi{j})
            Rtmp{VARi{j}(i)}.PolynomB(4) = 0;
        end
    end
end

for i=1:length(dipi)
    wdip(i) = Rtmp{dipi(i)}.PolynomB(2)/Rtmp{dipi(7)}.PolynomB(2);
end
for i=1:length(dipmi)
    wdipm(i) = Rtmp{dipmi(i)}.PolynomB(2)/Rtmp{dipmi(7)}.PolynomB(2);
end


if nargin < 1
    X0(1)  = atgetfieldvalues(Rtmp(qfii(1)),  'PolynomB',{1,2});
    X0(2)  = atgetfieldvalues(Rtmp(qfoi(1)),  'PolynomB',{1,2});
    X0(3)  = atgetfieldvalues(Rtmp(qfmi(1)),  'PolynomB',{1,2});
    X0(4)  = atgetfieldvalues(Rtmp(qfendi(1)),'PolynomB',{1,2});
    X0(5)  = atgetfieldvalues(Rtmp(qdendi(1)),'PolynomB',{1,2});
    X0(6)  = atgetfieldvalues(Rtmp(dipi(7)),  'PolynomB',{1,2});
    X0(7)  = atgetfieldvalues(Rtmp(dipmi(7)), 'PolynomB',{1,2});
    X0(8)  = atgetfieldvalues(Rtmp(sdi(1)),   'PolynomB',{1,3});
    X0(9)  = atgetfieldvalues(Rtmp(sdendi(1)),'PolynomB',{1,3});
    X0(10) = atgetfieldvalues(Rtmp(sfmi(1)),  'PolynomB',{1,3});
    X0(11) = atgetfieldvalues(Rtmp(sfoi(1)),  'PolynomB',{1,3});
    X0(12) = atgetfieldvalues(Rtmp(sfii(1)),  'PolynomB',{1,3});
end


VARi = {qfii; qfoi; qfmi; qfendi; qdendi};
for j = 1:5
    for i=1:length(VARi{j})
        Rtmp{VARi{j}(i)}.PolynomB(2) = X0(j);
        Rtmp{VARi{j}(i)}.K           = X0(j);
    end
end

clear VARi;
VARi = {dipi; dipmi};
W    = {wdip; wdipm};
for j = 1:2
    for i=1:length(VARi{j})
        Rtmp{VARi{j}(i)}.PolynomB(2) = X0(j+5) * W{j}(i);
        Rtmp{VARi{j}(i)}.K           = X0(j+5) * W{j}(i);
    end
end

clear VARi
VARi = {sdi; sdendi; sfmi; sfoi; sfii};
for j = 1:5
    for i=1:length(VARi{j})
        Rtmp{VARi{j}(i)}.PolynomB(3) = X0(j+7);
        Rtmp{VARi{j}(i)}.K           = X0(j+7);
    end
end

% ----------------------------
% force chromaticity to be 1/1
% ----------------------------
changechro=1;
if changechro == 1
    try
        disp('fitting chromaticity to [1,1]')
        Rtmp = atfitchrom(Rtmp,[1 1],'SFM','SDEND');
    catch
        disp('cannot fit chromaticity ... mmachine unstable')
    end
end

Rout = Rtmp;

end

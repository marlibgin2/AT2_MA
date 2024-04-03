function check_magnet_medmax(rrr, X)


qfii    = findcells(rrr,'FamName','QFI');
qfoi    = findcells(rrr,'FamName','QFO');
qfmi   = findcells(rrr,'FamName','QFM');
qfendi = findcells(rrr,'FamName','QFEND');
qdendi = findcells(rrr,'FamName','QDEND');

sdi     = findcells(rrr,'FamName','SD');
sdendi  = findcells(rrr,'FamName','SDEND');
sfmi    = findcells(rrr,'FamName','SFM');
sfoi    = findcells(rrr,'FamName','SFO');
sfii    = findcells(rrr,'FamName','SFI');

dipi   = findcells(rrr,'FamName','DIP');
dipmi  = findcells(rrr,'FamName','DIPm');
for i=1:length(dipi)
    wdip(i) = rrr{dipi(i)}.PolynomB(2)/rrr{dipi(7)}.PolynomB(2);
end
for i=1:length(dipmi)
    wdipm(i) = rrr{dipmi(i)}.PolynomB(2)/rrr{dipmi(7)}.PolynomB(2);
end


if ~isempty(X)
    VARi = {qfii; qfoi; qfmi; qfendi; qdendi};

    for j = 1:5
        for i=1:length(VARi{j})
            rrr{VARi{j}(i)}.PolynomB(2) = X(j);
            rrr{VARi{j}(i)}.K           = X(j);
        end
    end

    clear VARi;
    VARi = {dipi; dipmi};
    W    = {wdip; wdipm};
    for j = 1:2
        for i=1:length(VARi{j})
            rrr{VARi{j}(i)}.PolynomB(2) = X(j+5) * W{j}(i);
            rrr{VARi{j}(i)}.K           = X(j+5) * W{j}(i);
        end
    end

    VARi = {sdi; sdendi; sfmi; sfoi; sfii};
    for j = 1:5
        for i=1:length(VARi{j})
            rrr{VARi{j}(i)}.PolynomB(3) = X(j+7);
            rrr{VARi{j}(i)}.K           = X(j+7);
        end
    end

end

QFI_K = rrr{qfii(1)}.PolynomB(2); dQFI_K = 100*(QFI_K - 4.03) / 4.03; 
QFI_g  = 3/0.2998 * QFI_K;
QFI_length = rrr{qfii(1)}.Length;

QFO_K = rrr{qfoi(1)}.PolynomB(2); dQFO_K = 100*(QFO_K - 4.03) / 4.03; 
QFO_g  = 3/0.2998 * QFO_K;
QFO_length = rrr{qfoi(1)}.Length;

QFM_K = rrr{qfmi(1)}.PolynomB(2); dQFM_K = 100*(QFM_K - 3.774)/ 3.774; 
QFM_g  = 3/0.2998 * QFM_K;
QFM_length = rrr{qfmi(1)}.Length;

QFEND_K = rrr{qfendi(1)}.PolynomB(2); dQFEND_K = 100*(QFEND_K - 3.654)/ 3.654;  
QFEND_g  = 3/0.2998 * QFEND_K;
QFEND_length = rrr{qfendi(1)}.Length;

QDEND_K = rrr{qdendi(1)}.PolynomB(2);  dQDEND_K =100*(QDEND_K  + 2.504)/ -2.504; 
QDEND_g  = 3/0.2998 * QDEND_K;
QDEND_length = rrr{qdendi(1)}.Length;

DIP_K   = rrr{dipi(7)}.PolynomB(2); dDIP_K =100*(DIP_K  + 0.9)/ -0.9; 
DIP_g   = 3/0.2998 * DIP_K;

DIPM_K   = rrr{dipmi(7)}.PolynomB(2); 
DIPM_g   = 3/0.2998 * DIPM_K;

SD_K2 = rrr{sdi(1)}.PolynomB(3);
SD_g2 = 3/0.2998 * SD_K2;
SD_length = rrr{sdi(1)}.Length;

SDEND_K2 = rrr{sdendi(1)}.PolynomB(3);
SDEND_g2 = 3/0.2998 * SDEND_K2;
SDEND_length = rrr{sdendi(1)}.Length;

SFM_K2 = rrr{sfmi(1)}.PolynomB(3);
SFM_g2 = 3/0.2998 * SFM_K2;
SFM_length = rrr{sfmi(1)}.Length;

SFI_K2 = rrr{sfii(1)}.PolynomB(3);
SFI_g2 = 3/0.2998 * SFI_K2;
SFI_length = rrr{sfii(1)}.Length;

SFO_K2 = rrr{sfoi(1)}.PolynomB(3);
SFO_g2 = 3/0.2998 * SFO_K2;
SFO_length = rrr{sfoi(1)}.Length;



MagName = {'QFI'; 'QFO'; 'QFM'; 'QFEND'; 'QDEND'; 'DIP'; 'DIPM'; 'SD'; 'SDEND'; 'SFM'; 'SFI'; 'SFO'};
GraName = {'K (m-2)','B"/2 (T/m)'};
KMag    = [QFI_K; QFO_K; QFM_K; QFEND_K; QDEND_K; DIP_K; DIPM_K; 0; 0; 0; 0; 0];
dKMag   = [dQFI_K; dQFO_K; dQFM_K; dQFEND_K; dQDEND_K; DIP_K; DIPM_K; 0; 0; 0; 0; 0];
GMag    = [QFI_g; QFO_g; QFM_g; QFEND_g; QDEND_g; DIP_g; DIPM_g; 0; 0; 0; 0; 0];
K2Mag   = [0; 0; 0; 0; 0; 0; 0; SD_K2; SDEND_K2; SFM_K2; SFI_K2; SFO_K2];
G2Mag   = [0; 0; 0; 0; 0; 0; 0; SD_g2; SDEND_g2; SFM_g2; SFI_g2; SFO_g2];
%T = table(KMag,GMag,'RowNames',MagName)
T = table('RowNames',MagName);
T.('K (m⁻²)')     = num2str(KMag,4);
T.('dK (%)')     = num2str(dKMag,4);
T.('B´ (T/m)')  = num2str(GMag,4);
T.('K2 (m⁻³)')    = num2str(K2Mag,5); 
T.('B"/2 (T/m²)') = num2str(G2Mag,5);
disp(T); 

% https://ch.mathworks.com/matlabcentral/answers/254690-how-can-i-display-a-matlab-table-in-a-figure
% Get the table in string form.
TString = evalc('disp(T)')
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

end
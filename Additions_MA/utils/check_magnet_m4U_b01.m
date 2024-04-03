function check_magnet_m4U_b01(rrr, X)


q1i    = findcells(rrr,'FamName','Q1');
q2i    = findcells(rrr,'FamName','Q2');
q3i    = findcells(rrr,'FamName','Q3');
q4i    = findcells(rrr,'FamName','Q4');
q5i    = findcells(rrr,'FamName','Q5');
q6i    = findcells(rrr,'FamName','Q6');
r1i     = findcells(rrr,'FamName','R1');
r2i     = findcells(rrr,'FamName','R2');

s1i    = findcells(rrr,'FamName','S1');
s2i    = findcells(rrr,'FamName','S2');
s3i    = findcells(rrr,'FamName','S3');
s4i    = findcells(rrr,'FamName','S4');
s5i    = findcells(rrr,'FamName','S5');

o1i    = findcells(rrr,'FamName','O1');
o2i    = findcells(rrr,'FamName','O2');
o3i    = findcells(rrr,'FamName','O3');

%dmi   = findcells(rrr,'FamName','Dm');
d1i   = findcells(rrr,'FamName','D1');
d2i   = findcells(rrr,'FamName','D2');
d3i    = findcells(rrr,'FamName','D3');

%for i=1:length(dmi)
%    wm(i) = rrr{dmi(i)}.PolynomB(2)/rrr{dmi(7)}.PolynomB(2);
%end

for i=1:length(d1i)
    w1(i) = rrr{d1i(i)}.PolynomB(2)/rrr{d1i(7)}.PolynomB(2);
end

for i=1:length(d2i)
    w2(i) = rrr{d2i(i)}.PolynomB(2)/rrr{d2i(7)}.PolynomB(2);
end

for i=1:length(d3i)
    w3(i) = rrr{d3i(i)}.PolynomB(2)/rrr{d3i(7)}.PolynomB(2);
end


if ~isempty(X)
    VARi = {q1i; q2i; q3i; q4i; q5i; q6i; r1i; r2i};

    for j = 1:length(VARi)
        for i=1:length(VARi{j})
            rrr{VARi{j}(i)}.PolynomB(2) = X(j);
            rrr{VARi{j}(i)}.K           = X(j);
        end
    end

    clear VARi;
%    VARi = {dmi; d1i; d2i; d3i};
    VARi = {d1i; d2i; d3i};
    W    = {w1; w2; w3};
    for j = 1:length(VARi)
        for i=1:length(VARi{j})
            rrr{VARi{j}(i)}.PolynomB(2) = X(j+5) * W{j}(i);
            rrr{VARi{j}(i)}.K           = X(j+5) * W{j}(i);
        end
    end

    VARi = {s1i; s2i; s3i; s4i; s5i};
    for j = 1:length(VARi)
        for i=1:length(VARi{j})
            rrr{VARi{j}(i)}.PolynomB(3) = X(j+7);
            rrr{VARi{j}(i)}.K           = X(j+7);
        end
    end

end

Q1_K = rrr{q1i(1)}.PolynomB(2); %dQ1_K = 100*(Q1_K - 4.03) / 4.03; 
Q1_g  = 3/0.2998 * Q1_K;
Q1_length = rrr{q1i(1)}.Length;

Q2_K = rrr{q2i(1)}.PolynomB(2); dQ2_K = 100*(Q2_K - 4.03) / 4.03; 
Q2_g  = 3/0.2998 * Q2_K;
Q2_length = rrr{q2i(1)}.Length;

Q3_K = rrr{q3i(1)}.PolynomB(2); dQ3_K = 100*(Q3_K - 3.774)/ 3.774; 
Q3_g  = 3/0.2998 * Q3_K;
Q3_length = rrr{q3i(1)}.Length;

Q4_K = rrr{q4i(1)}.PolynomB(2); dQ4_K = 100*(Q4_K - 3.654)/ 3.654;  
Q4_g  = 3/0.2998 * Q4_K;
Q4_length = rrr{q4i(1)}.Length;

Q5_K = rrr{q5i(1)}.PolynomB(2); dQ5_K = 100*(Q5_K - 3.654)/ 3.654;  
Q5_g  = 3/0.2998 * Q5_K;
Q5_length = rrr{q5i(1)}.Length;

Q6_K = rrr{q6i(1)}.PolynomB(2); dQ6_K = 100*(Q6_K - 3.654)/ 3.654;  
Q6_g  = 3/0.2998 * Q6_K;
Q6_length = rrr{q6i(1)}.Length;

R1_K = rrr{r1i(1)}.PolynomB(2);  dR1_K =100*(R1_K  + 2.504)/ -2.504; 
R1_g  = 3/0.2998 * R1_K;
R1_length = rrr{r1i(1)}.Length;

R2_K = rrr{r2i(1)}.PolynomB(2);  dR2_K =100*(R2_K  + 2.504)/ -2.504; 
R2_g  = 3/0.2998 * R2_K;
R2_length = rrr{r2i(1)}.Length;

% DM_K   = rrr{dmi(7)}.PolynomB(2); dDM_K =100*(DM_K  + 0.9)/ -0.9; 
% DM_g   = 3/0.2998 * DM_K;

D1_K   = rrr{d1i(7)}.PolynomB(2); dD1_K =100*(D1_K  + 0.9)/ -0.9; 
D1_g   = 3/0.2998 * D1_K;

D2_K   = rrr{d2i(7)}.PolynomB(2); dD2_K =100*(D2_K  + 0.9)/ -0.9; 
D2_g   = 3/0.2998 * D2_K;

D3_K   = rrr{d3i(7)}.PolynomB(2); 
D3_g   = 3/0.2998 * D3_K;

S1_K2 = rrr{s1i(1)}.PolynomB(3);
S1_g2 = 3/0.2998 * S1_K2;
S1_length = rrr{s1i(1)}.Length;

S2_K2 = rrr{s2i(1)}.PolynomB(3);
S2_g2 = 3/0.2998 * S2_K2;
S2_length = rrr{s2i(1)}.Length;

S3_K2 = rrr{s3i(1)}.PolynomB(3);
S3_g2 = 3/0.2998 * S3_K2;
S3_length = rrr{s3i(1)}.Length;

S4_K2 = rrr{s4i(1)}.PolynomB(3);
S4_g2 = 3/0.2998 * S4_K2;
S4_length = rrr{s5i(1)}.Length;

S5_K2 = rrr{s5i(1)}.PolynomB(3);
S5_g2 = 3/0.2998 * S5_K2;
S5_length = rrr{s5i(1)}.Length;

O1_K3 = rrr{o1i(1)}.PolynomB(4);
O1_g3 = 3/0.2998 * O1_K3;
O1_length = rrr{o1i(1)}.Length;

O2_K3 = rrr{o2i(1)}.PolynomB(4);
O2_g3 = 3/0.2998 * O2_K3;
O2_length = rrr{o2i(1)}.Length;

O3_K3 = rrr{o3i(1)}.PolynomB(4);
O3_g3 = 3/0.2998 * O3_K3;
O3_length = rrr{o3i(1)}.Length;


% DM_theta = 0;
% for i=1:12
%     DM_theta = DM_theta + rrr{dmi(i)}.BendingAngle * 180/pi;
% end
D1_theta = 0;
for i=1:12
    D1_theta = D1_theta + rrr{d1i(i)}.BendingAngle * 180/pi;
end
D2_theta = 0;
for i=1:12
    D2_theta = D2_theta + rrr{d2i(i)}.BendingAngle * 180/pi;
end
D3_theta = 0;
for i=1:12
    D3_theta = D3_theta + rrr{d3i(i)}.BendingAngle * 180/pi;
end
R1_theta = 0;
for i=1:1
    R1_theta = R1_theta + rrr{r1i(i)}.BendingAngle* 180/pi;
end
R2_theta = 0;
for i=1:1
    R2_theta = R2_theta + rrr{r2i(i)}.BendingAngle* 180/pi;
end



MagName = {'Q1'; 'Q2'; 'Q3'; 'Q4'; 'Q5'; 'Q6'; 'R1'; 'R2'; 'D1'; 'D2'; 'D3'; 'S1'; 'S2'; 'S3'; 'S4'; 'S5'; 'O1'; 'O2'; 'O3'};
KMag    = [Q1_K; Q2_K; Q3_K; Q4_K; Q5_K; Q6_K; R1_K; R2_K; D1_K; D2_K; D3_K; 0; 0; 0; 0; 0; 0; 0; 0];
GMag    = [Q1_g; Q2_g; Q3_g; Q4_g; Q5_g; Q6_g; R1_g; R2_g; D1_g; D2_g; D3_g; 0; 0; 0; 0; 0; 0; 0; 0];
K2Mag   = [   0;    0;    0;    0;    0;    0;    0;    0;    0;    0;    0; S1_K2; S2_K2; S3_K2; S4_K2; S5_K2; 0; 0; 0];
G2Mag   = [   0;    0;    0;    0;    0;    0;    0;    0;    0;    0;    0; S1_g2; S2_g2; S3_g2; S4_g2; S5_g2; 0; 0; 0];
K3Mag   = [   0;    0;    0;    0;    0;    0;    0;    0;    0;    0;    0;     0;     0;     0;     0;     0; O1_K3; O2_K3; O3_K3];

thetaMag = [  0;    0;    0;    0;    0;    0; R1_theta; R2_theta; D1_theta; D2_theta; D3_theta; 0; 0; 0; 0; 0; 0; 0; 0];
%T = table(KMag,GMag,'RowNames',MagName)
T = table('RowNames',MagName);
T.('K (m⁻²)')     = num2str(KMag,4);
T.('B´ (T/m)')  = num2str(GMag,4);
T.('K2 (m⁻³)')    = num2str(K2Mag,5); 
T.('B"/2 (T/m²)') = num2str(G2Mag,5);
T.('K3 (m⁻4)')    = num2str(K3Mag,5); 
T.('theta (deg)') = num2str(thetaMag,5);
disp(T); 

% https://ch.mathworks.com/matlabcentral/answers/254690-how-can-i-display-a-matlab-table-in-a-figure
% Get the table in string form.
TString = evalc('disp(T)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
% % % annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
% % %     'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

end
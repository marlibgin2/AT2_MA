
% -------------------------------------------------------------------------
% test function to execute a coherent shift of all the slice elements in D1
% -------------------------------------------------------------------------
RING = m4U_240114_b01_02_03_02__grd;

% assign zero shifts to all elements, creating the T1, T2 fields
allelemi = 1: length(RING);
RING     = atsetshift(RING,allelemi,0,0);

d1i =  findcells(RING,'FamName','D1');
d2i =  findcells(RING,'FamName','D2');
d3i =  findcells(RING,'FamName','D3');
dmi =  findcells(RING,'FamName','Dm'); 
di  = [d1i d2i d3i dmi]; dcol = [100 100 100]/255; % color codes for plotting

q1i  = findcells(RING,'FamName','Q1');
q2i  = findcells(RING,'FamName','Q2');
q3i  = findcells(RING,'FamName','Q3');
q4i  = findcells(RING,'FamName','Q4');
q5i  = findcells(RING,'FamName','Q5');
r1i  = findcells(RING,'FamName','R1');
r2i  = findcells(RING,'FamName','R2');
qi = [q1i q2i q3i q4i q5i]; qcol=[255 100 0]/255;
ri = [r1i r2i]; rcol=[255 0 200]/255;


clear D1 D2 D3 DM Gi 
% ------------------------
% create groups of dipoles per achromat 
% ------------------------
for jA =1:1
    for i = 1:2
        DM{jA,i} = dmi(1+(i-1)*12:12+(i-1)*12);
        D1{jA,i} = d1i(1+(i-1)*12:12+(i-1)*12);
        D2{jA,i} = d2i(1+(i-1)*12:12+(i-1)*12);
    end
end
for jA =1:1
    D3{jA,1} = d3i(1+(1-1)*12:12+(1-1)*12);
end




%
% 
%


mGs = findcells(RING,'FamName','GRDs');
mGe = findcells(RING,'FamName','GRDe');

for jA = 1:1
    for iG = 1:7
        Gi{jA,iG} = [mGs(iG):mGe(iG)];
    end
end

%DX = 50e-6* randn(numel(Gi{1,1}),1);
%
% move the 7 girders in A1
%
clear DX DY
gdxdy = 100e-6;
qdxdy = 10e-6;
ddxdy = 10e-6;
sdxdy = 10e-6;
odxdy = 10e-6;

for iA = 1:1
    for ig =1:7
        gDX{iA,ig} =  gdxdy * ones(numel(Gi{iA,ig}),1);%*randn(1) 
        gDY{iA,ig} = -gdxdy * ones(numel(Gi{iA,ig}),1);%*randn(1)
    end
end

qDX = qdxdy*randn(numel(qi),1);
qDY = qdxdy*randn(numel(qi),1);
rDX = qdxdy*randn(numel(ri),1);
rDY = qdxdy*randn(numel(ri),1);

for iA=1:1
    for i = 1:7
        RING = ataddshift(RING, Gi{1,i}, gDX{iA,i}, gDY{iA,i});
    end
end
RING = ataddshift(RING, D1{1,1}, 10e-6,  21e-6);
RING = ataddshift(RING, D1{1,2}, 15e-6, -10e-6);
RING = ataddshift(RING, D2{1,1},-20e-6, -10e-6);
RING = ataddshift(RING, D2{1,2}, 15e-6, -30e-6);
% 
RING = ataddshift(RING, D3{1,1}, 30e-6, 60e-6);
% 
RING = ataddshift(RING, DM{1,1}, -30e-6, -60e-6);
RING = ataddshift(RING, DM{1,2}, -20e-6, -40e-6);
% 
% 
RING = ataddshift(RING, qi, qDX, qDY);

RING = ataddshift(RING, ri, rDX, rDY);

% -------------------
% plot the result ...
% -------------------
for i= 1:length(RING)
    dx(i) = -RING{i}.T1(1); dy(i) = -RING{i}.T1(3);
end
figure(2); clf; hold on
plot(spos,dx,'o')
plot(spos,dy,'o')
plot(spos(di),dx(di),'o','color',dcol,'MarkerFaceColor',dcol)
plot(spos(qi),dx(qi),'o','color',qcol,'MarkerFaceColor',qcol)
plot(spos(ri),dx(ri),'o','color',rcol,'MarkerFaceColor',rcol)
plot(spos(di),dy(di),'o','color',dcol,'MarkerFaceColor',dcol)
plot(spos(qi),dy(qi),'o','color',qcol,'MarkerFaceColor',qcol)
plot(spos(ri),dy(ri),'o','color',rcol,'MarkerFaceColor',rcol)



%stem(spos,dy)
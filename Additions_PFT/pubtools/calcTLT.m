function TousLT = CalcTLT(varargin)
% Wrapper for Touschek lifetime calculation
% Calculates, if necessary linear optics parameters
% expects a full ring as input, but may choose only
% one achromat for LMA calculations
%
[RING,LatticeOptData]= getargs(varargin,[],[]);
Ib = getoption(varargin,'Ib',0.5/176);
verbosef  = any(strcmpi(varargin,'verbose'));
integrationmethod = getoption(varargin,'integrationmethod','integral');
abstol = getoption(varargin, 'AbsTol', 1.0e-16);
reltol = getoption(varargin, 'Relol', 1.0e-16);
Nperiods = getoption(varargin,'Nperiods',20);
LMAperiods = getoption(varargin,'LMAperiods',1);
kcoupl = getoption(varargin,'kcoupl',0.025);


energy=RING{1}.Energy;
Periodicity=RING{1}.Periodicity;
circumference = findspos(RING,length(RING)+1)*Periodicity;
period = circumference/Nperiods;
ats = atsummary(RING);
emitx = ats.naturalEmittance/(1+kcoupl);
emity = ats.naturalEmittance*kcoupl/(1+kcoupl);
sigp  = ats.naturalEnergySpread;
sigs  = ats.bunchlength;
deltalimit = ats.energyacceptance;

S0min = 0.0;
S0max = LMAperiods*period;
if (verbosef)
    [map_l,map_h,~,Ipos,~] = CalcPlotMA(RING,LatticeOptData,'lmafams','All','S0max',...
                             S0max,'S0min',S0min,'deltalimit',deltalimit,'verbose');
else
    [map_l,map_h,~,Ipos,~] = CalcPlotMA(RING,LatticeOptData,'lmafams','All','S0max',...
                             S0max,'S0min',S0min,'deltalimit',deltalimit);
end
lmap=cat(2,map_h,map_l);
[~,lindata] = atlinopt4(RING,Ipos);
for i=1:length(lindata)
    lindata(i).Length = RING{Ipos(i)}.Length;
end

if (verbosef)
    TousLT = TLT(lindata,lmap,'Ib',Ib,'circumference',circumference,'energy',energy,'emitx',emitx,'emity',emity,...
          'sig',sigp,'sigs',sigs,'abstol',abstol,'reltol', reltol, 'integrationmethod', integrationmethod,'verbose');
else
    TousLT = TLT(lindata,lmap,'Ib',Ib,'circumference',circumference,'energy',energy,'emitx',emitx,'emity',emity,...
          'sig',sigp,'sigs',sigs,'abstol',abstol,'reltol', reltol, 'integrationmethod', integrationmethod);
end


function plotPareto(MOGAResults)
%Plots ParetoFron from MOGA Results
% inout is MOGAREsults structure genrated by funciton MOGA.m
%
LatticeOptData = MOGAResults.LatticeOptData;
nvars          = LatticeOptData.nvars;
Trb            = LatticeOptData.Trb;
fitvars        = MOGAResults.fitvars;


if (isfield(LatticeOptData,'optMode'))
    optMode  = LatticeOptData.optMode;
else
    optMode = 'SIMP';
end

RunNumber = MOGAResults.RunNumber;
x=MOGAResults.ParetoFront(:,nvars+1);

%if strcmp(optMode,'CHRO')
%    y=1./MOGAResults.ParetoFront(:,nvars+2);
%else
%    y=MOGAResults.ParetoFront(:,nvars+2);
%end
y=MOGAResults.ParetoFront(:,nvars+2);
figure;plot(x,y,'o');
xlabel(fitvars{1});
if strcmp(optMode,'CHRO')
    
    ylabel('-AD [mm*2]');
else
    ylabel(fitvars{2});
end

title(strcat(RunNumber,'# ',optMode, '# RB= ',num2str(Trb)));
grid;
end


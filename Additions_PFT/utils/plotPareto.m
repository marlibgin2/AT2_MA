function plotPareto(varargin)
% Plots ParetoFront from MOGA Results
% input is a cell array of MOGAResults structures generated 
% by the function MOGA.m
%
MOGAResults    = getargs(varargin,[]);
legend_str     = getoption(varargin,'legend', '');

if (iscell(MOGAResults))
    nFronts = numel(MOGAResults);
    MOGAR   = MOGAResults{1};
else
    MOGAR   = MOGAResults; 
    nFronts = 1;
end

LatticeOptData = MOGAR.LatticeOptData;
nvars          = LatticeOptData.nvars;
fitvars        = MOGAR.fitvars;


if (isfield(LatticeOptData,'optMode'))
    optMode  = LatticeOptData.optMode;
else
    optMode = 'SIMP';
end

RunNumber{1} = MOGAR.RunNumber;
x=MOGAR.ParetoFront(:,nvars+1);
y=MOGAR.ParetoFront(:,nvars+2);

figure;plot(x,y,'o');
xlabel(fitvars{1});

if strcmp(optMode(1),'CHRO')
    ylabel('-AD [mm*2]');
else
    ylabel(fitvars{2});
end

title(strcat({'Optmization Mode = '}, optMode));
grid;
if (isempty(legend_str))
    legendstr{1}=strrep(RunNumber{1},'_','-');
else
    legendstr = legend_str;
end

if (nFronts>1)
    hold on;
    for i=2:nFronts
        nvars = MOGAResults{i}.LatticeOptData.nvars;
        if (isempty(legend_str))
            legendstr{i} = strrep(MOGAResults{i}.RunNumber,'_','-');
        end
        x=MOGAResults{i}.ParetoFront(:,nvars+1);
        y=MOGAResults{i}.ParetoFront(:,nvars+2);
        plot(x,y,'o');
    end
    
end

legend(legendstr);

function [X, fval, MOGAResults] = recoverMOGAdata(filin, X0, dx)
path = pwd; 
if nargin<1
    [filin, path] = uigetfile('nMOGAgen_*.mat');
end
if nargin<2
    dx = [0.075,     0.075,    0.075,     0.075,     0.075,     0.075,     0.075,    0.6*4,  0.6*4,  0.6*4,  0.6*4,  0.6*4];
    X0 = [5.3766    5.5046    4.2275    4.1357   -2.6794   -1.0479   -0.7410 -243.4325  -22.9891  166.8030  317.8298  296.1150];
end

date = filin(10:17);
time = filin(19:26);

FF = load([path '/' filin]);
goodrank = find(FF.state.Rank==1);

X    =  FF.state.Population(goodrank,:);
fval =  FF.state.Score(goodrank,:);

MOGAResults.options     = FF.options;
MOGAResults.FuncTol     = FF.options.FunctionTolerance;  
MOGAResults.ConstTol    = FF.options.ConstraintTolerance; 
MOGAResults.MaxGen      = FF.options.MaxGenerations;   
MOGAResults.output      = [];
MOGAResults.sizPop      = FF.options.PopulationSize;   
MOGAResults.ParetoFront = cat(2,X,fval); 
MOGAResults.X           = X;        
MOGAResults.fval        = fval;     
MOGAResults.exitflag    = FF.state.StopFlag; 
MOGAResults.lb          = FF.options.InitialPopulationRange(1,:);
MOGAResults.ub          = FF.options.InitialPopulationRange(2,:);
MOGAResults.X0          = X0;
MOGAResults.dx          = dx;
MOGAResults.population  = FF.state.Population;
MOGAResults.score       = FF.state.Score;
MOGAResults.nvars       = size(FF.state.Population,2);


output   = MOGAResults.output;
exitflag = 1; 
save(['MOGArecover_' date '_' time '.mat'],'X','fval','exitflag','output','MOGAResults')
figure(1010)
p = plot(fval(:,1), fval(:,2),'o');


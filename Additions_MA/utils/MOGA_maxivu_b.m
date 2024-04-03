% --------------------------------------------------------
% MOGA optimisation 
% ---------------------------------------------------------
function [X,fval,exitflag,output,MOGAResults] = MOGA_maxivu_b(X0)%, besti)

rrr = PrepareInitialRing(X0); % create a ring 'rrr' of type 'C' with quad/sext/oct as from X0

[MOGAfilename MOGAdir]=uigetfile('MOGA_*.mat');
if (ischar(MOGAfilename))
        load([MOGAdir MOGAfilename],'MOGAResults');
    else
        fprintf('No file chosen. Aborting... \n');
        return
    end

%% trying to buzz around 155 pm.rad / I reduce the quad span and increase the sextupole changes
%% delta  = [0.44/1,  0.44/1, 0.38/1,  0.38/1,  0.25/1,  0.09/1,  0.08/1, 11*6, 16*6, 17*6, 17*6, 20*6]*1; % 5% variation for quad/Dipgrad, 60% variation for sext
dx     = [0.09,     0.09,    0.09,     0.09,     0.09,     0.09,     0.09,    0.6*4,  0.6*4,  0.6*4,  0.6*4,  0.6*4];
nvars  = length(dx);
sizPop =  500; 
MaxGen   = 20;
FuncTol  = 1e-3;
ConstTol = 5e-1;

if (size(X0,2)==nvars)
    lb = X0-abs(X0).*dx/5000;
    ub = X0+abs(X0).*dx/5000;
    fprintf('lower bounds reset to \r\n');
    fprintf('upper bounds reset to \r\n');
end

if 1==0 %old way
    iniPop = 2*(rand(sizPop,12)-0.5).*delta/1+X0;
else % a la Pedro
    OldPopulation     = MOGAResults.population;
    OldPopSize        = MOGAResults.sizPop;
    Oldnvars          = MOGAResults.nvars;
    mvars             = min(nvars,Oldnvars);
    mPopSize          = min(OldPopSize,sizPop);
    iniPop = zeros(sizPop,nvars);
    for i=1:sizPop
          for j=1:nvars
               IniPop(i,j) = (lb(j)+ub(j))/2;
          end
    end
    IniPop(1:mPopSize,1:mvars)=OldPopulation(1:mPopSize,1:mvars);

end

% iniPop = [iniPop; X0]; 

%sizIni = 100;
%iniPop = 2*(rand(sizIni,12)-0.5).*delta/100+X0;
%iniPop =  repmat(X0,sizPop,1);


% iniPop = 2*(rand(5000-size(X0,1),12)-0.5).*delta+X0(besti,:);
% iniPop = [iniPop; X0]; 

%iniPop = X0; 
controlgcp = 1; 
usepar = true;
if controlgcp==1
   c = parcluster;
   % next line is optional at MAX IV
   c.AdditionalProperties.AccountName = 'any-virtual-account-name';
   % 6 hour walltime
   c.AdditionalProperties.WallTime = '70:00:00';
   % hyperthreading enabled
   c.saveProfile
   c.NumThreads = 2;
   parpool('aurora R2022a',128) %48 %96
% % % %    c.NumThreads = 1;
% % % %    parpool('local',12) %48
   pp = gcp; 
   usepar = true;
end



%
% start with the initial options
%
options = optimoptions('gamultiobj');

options = optimoptions(options,'PlotFcn',@gaplotpareto);
%options = optimoptions(options,'display', 'iter');
options = optimoptions(options,'PopulationSize',sizPop);
options = optimoptions(options,'InitialPopulationMatrix',iniPop);
%options = optimoptions(options,'InitialPopulation',iniPop);
options = optimoptions(options,'UseParallel',usepar);
options = optimoptions(options,'UseVectorized', false);
options = optimoptions(options,'OutputFcns', @MOGAoutputfcn);
options = optimoptions(options,'FunctionTolerance',FuncTol);%1e-3);
options = optimoptions(options,'MaxGenerations',MaxGen);%40
options = optimoptions(options,'ConstraintTolerance',ConstTol);%1e-1);


% options = optimoptions('gamultiobj','PlotFcn',@gaplotpareto,'PopulationSize',sizPop,'InitialPopulationMatrix',iniPop,...
%                        'UseParallel',usepar,'UseVectorized', false,...
%                        'display', 'iter', 'OutputFcns', @MOGAoutputfcn, ...
%                           'FunctionTolerance',1e-3,'MaxStallGenerations',5,'MaxGenerations',20, ...
%                        'ConstraintTolerance',1e-1);
%%% [X,fval,exitflag,output] = fminsearchbnd(@fun_medmax_match_AT2,X0, X0-delta, X0+delta, options);X
%%% [X,fval,exitflag,output] = fmincon(@fun_medmax_match_AT2,X0, [], [], [], [], X0-delta, X0+delta, [], options);

%fitnessfcn = @(X)[obj_emix(X, rrr), obj_rdt(X, rrr) ];
%%%%% fitnessfcn = @(X)[obj_emix(X, rrr), obj_modelDA(X, rrr) ];
fitnessfcn = @(X)OBJfcn(X, rrr); % object function
mycon      = @(X)CONfcn(X, rrr); % constraint function
nvars      = length(X0);
%[X,fval,exitflag,output,population,score] = gamultiobj(fitnessfcn,nvars,[],[],[],[],X0-delta, X0+delta, mycon, options);
[X,fval,exitflag,output,population,score] = gamultiobj(fitnessfcn,nvars,[],[],[],[],lb, ub, mycon, options);

 MOGAResults.options     = options;
 MOGAResults.FuncTol     = FuncTol;  %
 MOGAResults.ConstTol    = ConstTol; %
 MOGAResults.MaxGen      = MaxGen;   %
 MOGAResults.output      = output;
 MOGAResults.sizPop      = sizPop;   %
 MOGAResults.ParetoFront = cat(2,X,fval); %
 MOGAResults.X           = X;        %
 MOGAResults.fval        = fval;     %
 MOGAResults.exitflag    = exitflag; %
 MOGAResults.lb          = lb;
 MOGAResults.ub          = ub;
 MOGAResults.X0          = X0;
 MOGAResults.dx          = dx;
 MOGAResults.population  = population;
 MOGAResults.score       = score;
 MOGAResults.nvars       = nvars;
adesso   = datestr(now,'ddmmyyyy_HHMMSS');
filename = ['MOGA_' adesso '.mat'];
%MOGAResults.RunNumber   = filename;

save(filename,'X','fval','exitflag','output','MOGAResults')
%save(filename,'MOGAResults')


if controlgcp==1
    delete(pp)
end

end


% --------------------------------------------------------
% MOGA optimisation of emittance and dynamic aperture
% ---------------------------------------------------------
function [X,fval,exitflag,output,MOGAResults] = MOGA_TEST(X0)%, besti)
if nargin<1
    X0 = [];
end

dx     = [0.10, 0.10, 0.10, 0.10, 0.10, 0.10, ... % QUAD
          0.10, 0.10, ...                         % RB quad   
          2.00, 2.00, ...                         % BB quad 
          2.00, 2.00, 2.00, 2.00, 2.00, ...       % SEXT   
          2.00, 2.00, 2.00]*1e-3;                % OCT  
nvars  = length(dx);
sizPop   = 5000;  % 5000
MaxGen   = 100; % 20
FuncTol  = 1e-3;
ConstTol = 5e-1;

if (size(X0,2)==nvars)
    lb = X0-abs(X0).*dx/0.001;
    ub = X0+abs(X0).*dx/10000;
    fprintf('lower bounds reset to \r\n');
    fprintf('upper bounds reset to \r\n');
end

if 1==1 %old way
    iniPop = 2*(rand(sizPop,nvars)-0.5).*dx+X0;
else % a la Pedro
    [MOGAfilename MOGAdir]=uigetfile('../MOGA_Data/MOGA_*.mat');
    if (ischar(MOGAfilename))
        load([MOGAdir MOGAfilename],'MOGAResults');
    else
        fprintf('No file chosen. Aborting... \n');
        return
    end
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
   parpool('aurora R2022a',256); % 128 %48 %96
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
options = optimoptions(options,'PlotInterval',1);
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

myobj      = @(X)OBJfcn_TEST(X); % object function
%mycon      = @(X)CONfcn_TEST(X); % constraint function
nvars      = length(X0);
%[X,fval,exitflag,output,population,score] = gamultiobj(myobj,nvars,[],[],[],[],lb, ub, mycon, options);
[X,fval,exitflag,output,population,score] = gamultiobj(myobj,nvars,[],[],[],[],lb, ub, options);

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
filename = ['TEST_MOGA_' adesso '.mat'];
%MOGAResults.RunNumber   = filename;

save(filename,'X','fval','exitflag','output','MOGAResults')
%save(filename,'MOGAResults')


if controlgcp==1
    delete(pp)
end

end


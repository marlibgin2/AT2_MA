function MOGAResults = MOGA(varargin)
% Multiobjective optimization of lattice parameters using genetic
% algorithm
%
%% Inputs
% Mandatory arguments:
%         LatticeOptData : structure withOptimzation configiurationinfo
%         optfnct : function to be minimized
%              LattOpt_EmitCHRO  : Emittance and squared sum of chromaticities 
%              LattOpt_EmitDynap : Emittance and on-momentum Dynamic Aperture
%              LattOpt_EmitRDT   : Emittance and Resonant Driving Terms
%
% Optional arguments input with syntax ('ParameterName', parameter value)
%
%         X0: initial guess - if = [], the lower an upper bounds given by
%             the lb and ub input parameters are passed on to the genetic
%             optmizer. Otherwise, the lb and ub are calculated from the
%             initial guess +/- dx*initialguess/100.
%         dx: relative change of decision variables to define the lb and ub
%             parameters. Ony relevant if X0 is not []
%         lb: lower bounds (1xnvars) array, where nvars is the number of variable to vary
%         ub: upper bounds (1xnvars) array, where nvars is the number of variable to vary
%         PopSize: population size
%         FuncTol: Tolerance for stopping MOGA generations
%         MaxGens: Maximum number of generations
%         X0_Lin : Linear lattice decision variables - only relevant for
%                  "SEXT" (legacy) or "NonLinear" optmization mode.
% 
% Optional flags (if present the option is activated)
%
%         save  : if present saves out structure onto a file in folder
%                "\MOGA_Scans" and automatically generated file name.
%         cont  : if present assumes this is a continuation run and uses the
%                final population of a previus run (selected from a dialog box) as
%                the initial population for this run. Otherwise the initial population is generated acording
%                to the choice indicated below in "options"
%         comp  : if present, collects additional ParetoFront data
%
%% Outputs
% MOGAResults structure with fields:
%
%
%% Usage examples
% MOGAResults=MOGA(LatticeOptData,'RINGOpt_EmitDynAp','X0',X0_b3,'dx',20, 'PopSize',1000,'MaxGens',50,'comp','save','cont');
% MOGAResults=MOGA(LatticeOptData,'RINGOpt_EmitDynAp','X0',X0_b3,'dx',20, 'PopSize',1000,'MaxGens',50,'comp');

%% Input Argument Parsing
[LatticeOptData, optfnct]=getargs(varargin,[],'LattOpt_EmitDynAp');
X0             = getoption(varargin,'X0',[]);
dx             = getoption(varargin,'dx',[]);
PopSize        = getoption(varargin,'PopSize',1000);
FuncTol        = getoption(varargin,'FuncTol', 1E-9);
MaxGens        = getoption(varargin,'MaxGens', 100);
X0_Lin         = getoption(varargin,'X0_Lin', []);
lb             = getoption(varargin,'lb',[]);
ub             = getoption(varargin,'ub',[]);

savef          = any(strcmpi(varargin,'save'));
contf          = any(strcmpi(varargin,'cont'));
compf          = any(strcmpi(varargin,'comp'));

%% Preamble
%
nvars     = LatticeOptData.nvars;
scan_fams = LatticeOptData.scan_fams; 
lattMode  = LatticeOptData.lattMode;
optMode   = LatticeOptData.optMode;

switch optfnct
    case 'LattOpt_EmitChro'
        fitvars     = {'Emit[pmrad]';'Sqrt(Chrox**2+Chroy**2)'};
    case {'LattOpt_EmitDynAp';'RINGOpt_EmitDynAp'}
        fitvars     = {'Emit[pmrad]';'-AD[mm**2]'};  
    case 'LattOpt_EmitRDT'
        fitvars     = {'Emit[pmrad]';'RDT'};
    case 'LattOpt_EmitDiffRate'
        fitvars     = {'Emit[pmrad]';'DiffusionRate'};

    otherwise    
        fprintf('Invalid optimization function %s . Aborting .\n', optfnct);
        MOGAResults=[];
        return;
end
fh = str2func(optfnct);
LattOptfnct = @(x)fh(x,LatticeOptData);
%
% Checks compatibility between optimization mode and optimization function
%
switch optMode
    case {'Linear'}
       if (not(strcmp(optfnct,'LattOpt_EmitChro')))
           fprintf('Incompatible optimization function %s \n', optfnct);
           fprintf('for optimization mode %s \n', optMode);
           fprintf('Aborting...\n');
           MOGAResults=[];
           return
       end
    case {'NonLinear','Full','FullOct'}
      if (strcmp(optfnct,'LattOpt_EmitChro'))
           fprintf('Incompatible optimization function %s \n', optfnct);
           fprintf('for optimization mode %s \n', optMode);
           fprintf('Aborting...\n');
           MOGAResults=[];
           return
      end
   otherwise
       fprintf ('Invalid optimization mode. Aborting \n');
       MOGAResults=[];
       return
end
% 
%
fprintf('----- \n');
fprintf('%s Starting MOGA Scan \n', datetime);
fprintf('Lattice Mode      = %s  \n', lattMode);
fprintf('Optimization Mode      = %s  \n', optMode);
fprintf('Optimization Function  = %s  \n', optfnct);
fprintf('Decision Variables     = %s ,%s ,%s ,%s ,%s ,%s ,%s ,%s ,%s ,%s \n', ...
    scan_fams{1:length(scan_fams)});
fprintf(' \n ----- \n');

if (size(X0,2)==nvars)
    lb = X0-abs(X0)*dx/100;
    ub = X0+abs(X0)*dx/100;
end

fprintf('lower bounds set to: ');
fprintf('%5.3f ', lb(1:nvars));
fprintf('\r\n');
fprintf('upper bounds set to: ');
fprintf('%5.3f ', ub(1:nvars));
fprintf('\r\n');

% Start with the default options
options = optimoptions('gamultiobj');

% Modify options setting

% creation options
%Range = cat(1,lb,ub);
%options = optimoptions(options,'CreationFcn', @gacreationuniform);
%options = optimoptions(options,'InitialPopulationRange', Range);
options = optimoptions(options,'CreationFcn', @gacreationnonlinearfeasible);
options = optimoptions(options,'MutationFcn', @mutationadaptfeasible);
options = optimoptions(options,'CrossoverFcn',@crossoverintermediate);
options = optimoptions(options,'Display', 'iter');
options = optimoptions(options,'FunctionTolerance', FuncTol);
options = optimoptions(options,'MaxGenerations', MaxGens);
options = optimoptions(options,'PopulationSize', PopSize);
options = optimoptions(options,'UseParallel', true, 'UseVectorized', false);
%options = optimoptions(options,'PlotFcn', { @gaplotpareto, @gaplotrankhist });
options = optimoptions(options,'PlotFcn', { @gaplotpareto });

% if this is a continuation run, get the initial population from a previous
% MOGA result

if(contf)
    filename=uigetfile;
    if (ischar(filename))
        load(filename,'MOGAResults');
        if(isfield(MOGAResults.LatticeOptData,'optMode'))
            OldoptMode=MOGAResults.LatticeOptData.optMode;
        else
            OldoptMode='SIMP';
        end
    else
        fprintf('No file chosen. Aborting... \n');
        return
    end
  
    switch optMode
        case 'Linear' 
            cancontinue = ismember(OldoptMode,['Linear']);

        case 'NonLinear'
            cancontinue = ismember(OldoptMode,['NonLinear']);

        case 'Full'
            cancontinue = ismember(OldoptMode,['Full']);

        case 'FullOct'
            cancontinue = ismember(OldoptMode,['Full' 'FullOct']);
    end
    if(not(cancontinue)) 
        fprintf('Inconsistent optmization modes (%s to %s ) in continuation run. Aborting...\n',...
              OldoptMode, optMode);
        return
    else
        OldPopulation     = MOGAResults.population;
        OldPopSize        = MOGAResults.PopSize;
        Oldnvars          = MOGAResults.LatticeOptData.nvars;
        mvars             = min(nvars,Oldnvars);
        mPopSize          = min(OldPopSize,PopSize);
        InitialPopulation = zeros(PopSize,nvars);
        for i=1:PopSize
            for j=1:nvars
                InitialPopulation(i,j) = (lb(j)+ub(j))/2;
            end
        end
        InitialPopulation(1:mPopSize,1:mvars)=OldPopulation(1:mPopSize,1:mvars);
        if(contf&&isfield(MOGAResults,'DVLins'))
            DVLins=MOGAResults.DVLins;
        else
            DVLins = X0_Lin;
        end
    end 
    options = optimoptions(options,'InitialPopulation', InitialPopulation);
    MOGAResults=[];
    MOGAResults.ContinuationRun='Yes';
    MOGAResults.PreviousRun=filename;
    fprintf ('Continuation run from:  %s \n', filename);
    fprintf ('Optmization Modes:      %s to %s \n', OldoptMode, optMode);
    fprintf ('Population:             %4d to %4d \n', OldPopSize, PopSize);
    fprintf ('Decision Variables:     %4d to %4d \n', Oldnvars, nvars);
else
    MOGAResults.ContinuationRun='No';
    MOGAResults.PreviousRun='';
    DVLins = X0_Lin;
end

% sets fixed linear lattice parameters in case of pure Nonlinear
% optimization case 
%
if (strcmp(optMode,'NonLinear'))
  nlinfams =  size(LatticeOptData.lin_fams,1);
  if (nlinfams~=size(DVLins,2))
    fprintf('Linear Optics Vector size = %2d does not match lin fams vector size = %d , aborting... \n', size(DVLins,2), nlinfams)
    return
  else
    LatticeOptData.HACHRO  = setLins(1,LatticeOptData.HACHRO,LatticeOptData,DVLins);
    LatticeOptData.ACHRO   = setLins(2,LatticeOptData.ACHRO,LatticeOptData,DVLins);
    LatticeOptData.UC      = setLins(3,LatticeOptData.UC,LatticeOptData,DVLins);
    LatticeOptData.RING    = setLins(5,LatticeOptData.RING,LatticeOptData,DVLins);
    LatticeOptData.RINGGRD = setLins(6,LatticeOptData.RINGGRD,LatticeOptData,DVLins);
  end
end

Constfnct=@(x)stabcon(x,LatticeOptData);

%% Runs genetic algorithm
tstart=tic;

try
    [x,fval,exitflag,output,population,score] = gamultiobj...
    (LattOptfnct,nvars,[],[],[],[],lb,ub,Constfnct,options);
catch ME
    fprintf('%s Error during MOGA run \n', datetime);
    fprintf('Error message was:%s \n',ME.message);
end

%% Collects results in output structure
MOGAResults.LatticeOptData=LatticeOptData;
MOGAResults.options=options;
MOGAResults.opffnct=optfnct;
MOGAResults.fitvars=fitvars;
MOGAResults.constraints = {'abs(Tr(Mx))<2';'abs(Tr(My))<2';...
                    strcat('BetaX<',num2str(LatticeOptData.BetaXMAX));...
                    strcat('BetaY<',num2str(LatticeOptData.BetaYMAX));...
                    strcat('EtaX0<',num2str(LatticeOptData.EtaX0MAX));...
                    strcat('abs(FractuneGoal_X-FracTuneX) < ',num2str(LatticeOptData.DeltaFracTuneX));...
                    strcat('abs(FractuneGoal_Y-FracTuneY) < ',num2str(LatticeOptData.DeltaFracTuneY));...
                    strcat('BetaY0 - BetaX0 < ', num2str(LatticeOptData.DBetaXY));...
                    strcat('BetaX0 < ', num2str(LatticeOptData.BetaX0Max));...
                    strcat('BetaY0 < ', num2str(LatticeOptData.BetaY0Max));...
                    strcat('Jx < 3.0 ')};
MOGAResults.FuncTol=FuncTol;
MOGAResults.MaxGens=MaxGens;
MOGAResults.PopSize=PopSize;
MOGAResults.ParetoFront=cat(2,x,fval);
MOGAResults.exitflag=exitflag;
MOGAResults.output=output;
MOGAResults.population=population;
MOGAResults.score=score;
MOGAResults.lb=lb;
MOGAResults.ub=ub;
MOGAResults.X0=X0;
MOGAResults.dx=dx;
filename = strcat('MOGA_',datestr(now,30));
MOGAResults.RunNumber = filename;
if(strcmp(optMode,'NonLinear'))
    MOGAResults.DVLins=DVLins;
end

telapsed = toc(tstart);
MOGAResults.telapsed = telapsed;


nfitvars = size(fitvars,1);
%
% complements lattice data for Pareto Front Indviduals
%
if (compf)
   fprintf('Colllecting additional ParetoFront data \n'); 
   pfsize=size(MOGAResults.ParetoFront,1); 
   tic;
   for i=1:pfsize
       try
        if (strcmpi(optMode,'NonLinear')||strcmpi(optMode,'Full')||strcmpi(optMode,'FullOct'))
            rpar=ExMOGA(MOGAResults,-i,'fitchrom');
        else
            rpar=ExMOGA(MOGAResults,-i);
        end 
        MOGAResults.ParetoFront(i,nvars+nfitvars+1)=rpar.outputs.tunesuc0(1);
        MOGAResults.ParetoFront(i,nvars+nfitvars+2)=rpar.outputs.tunesuc0(2);
        MOGAResults.ParetoFront(i,nvars+nfitvars+3)=rpar.outputs.atsummary.beta0(1);
        MOGAResults.ParetoFront(i,nvars+nfitvars+4)=rpar.outputs.atsummary.beta0(2);
        MOGAResults.ParetoFront(i,nvars+nfitvars+5)=rpar.outputs.atsummary.etax;
        MOGAResults.ParetoFront(i,nvars+nfitvars+6)=rpar.outputs.atsummary.damping(1);
        MOGAResults.ParetoFront(i,nvars+nfitvars+7)=rpar.outputs.Sc1;
        MOGAResults.ParetoFront(i,nvars+nfitvars+8)=rpar.outputs.Sc2;
        MOGAResults.ParetoFront(i,nvars+nfitvars+9)=i;
       catch ME
          fprintf('Error getting linear optics for individual n. %4d \n', i);
          fprintf('Error message was:%s \n',ME.message);
       end
   end
   toc;
end
%
%% Saves results file

if(savef)
    fprintf('Saving file %s \n', filename);
    filename=strcat('/home/pedtav/Documents/Max4U/MOGA_Scans/',filename);
    try
        save(filename,'MOGAResults');
    catch
        fprintf('%s Problems saving MOGA Results file', datetime)
    end
end

fprintf('----- \n');
fprintf('MOGA Scan ended on %s  \n', datetime);
fprintf('----- \n');



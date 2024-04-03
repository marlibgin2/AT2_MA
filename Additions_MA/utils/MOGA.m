function MOGAResults = MOGA(optfnct,X0, dx, lb,ub, PopSize, FuncTol, MaxGens, savef, contf, compf, LatticeOptData, X0_Lin)
% Runs Multiobjective optimization of lattice parameters using genetic
% algorithm
%
% inputs:
%         optfnct : function to be minimized
%              LattOpt_EmitCHRO  : Emittance and squared sum of chromaticities 
%              LattOpt_EmitDynap : Emittance and on-momentum Dynamic Aperture
%              LattOpt_EmitRDT   : Emittance and Resonant Driving Terms
%
%         X0: initial guess - if = [], the lower an upper bounds given by
%             the lb and ub input parameters are passed on to the genetic
%             optmizer. Otherwisee, the lb and ub are calculated from the
%             initial guess +/- dx*initialguess/100.
%         dx: relative change of decision variables to define the lb and ub
%             parameters. Ony relevant if X0 is not []
%         lb: lower bounds (1xnvars) array, where nvars is the number of variable to vary
%         ub: upper bounds (1xnvars) array, where nvars is the number of variable to vary
%         PopSize: population size
%         FuncTol: Tolerance for stopping MOGA generations
%         MAXGens: Maximum number of generations
%         savef  : if 'Y', saves out structure onto a file in folder
%                "\MOGA_Scans" and automatically generated file name.
%         contf  : if 'Y' assumes this is a continuation run and uses the
%                final population of a previus run (selected from a dialog box) as
%                the initial population for this run. Otherwise the initial population is generated acording
%                to the choice indicated below in "options"
%         compf  : collects additional ParetoFront data
%         X0_Lin : Linear lattice decision variables - only relevant for
%                  "SEXT" optmization mode.
% 
%% NOTE 1: Various options of the minimization algorithm can be chosen 
%% by modifying the calls to optimoptions below
%
%

nvars     = LatticeOptData.nvars;
scan_fams = LatticeOptData.scan_fams; 
optMode   = LatticeOptData.optMode;
revBmode  = LatticeOptData.revBmode;


switch optfnct
    case 'LattOpt_EmitChro'
        fitvars     = {'Emit[pmrad]';'Sqrt(Chrox**2+Chroy**2)'};
    case 'LattOpt_EmitDynAp'
        fitvars     = {'Emit[pmrad]';'-AD[mm**2]'};  
    case 'LattOpt_EmitRDT'
        fitvars     = {'Emit[pmrad]';'RDT'};
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
   case {'SIMP','COMP'}
       if (not(strcmp(optfnct,'LattOpt_EmitChro')))
           fprintf('Incompatible optimization function %s \n', optfnct);
           fprintf('for optimization mode %s \n', optMode);
           fprintf('Aborting...\n');
           MOGAResults=[];
           return
       end
   case {'CHRO','SEXT'}
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

fprintf('----- \n');
fprintf('%s Starting MOGA Scan \n', datetime);
fprintf('Optimization Mode      = %s  \n', optMode);
fprintf('Reverse Bend Mode      = %s  \n', revBmode);
fprintf('Optimization Function  = %s  \n', optfnct);
fprintf('Decision Variables     = %s ,%s ,%s ,%s ,%s ,%s ,%s ,%s ,%s ,%s \n', ...
    scan_fams{1:length(scan_fams)});
fprintf(' \n ----- \n');

if (size(X0,2)==nvars)
    lb = X0-abs(X0)*dx/100;
    ub = X0+abs(X0)*dx/100;
    fprintf('lower bounds reset to \r\n');lb
    fprintf('upper bounds reset to \r\n');ub
end

%% Start with the default options
options = optimoptions('gamultiobj');

%% Modify options setting

% creation options
Range = cat(1,lb,ub);
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

if(strcmp(contf,'Y'))
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
  
    if( (strcmp(OldoptMode,'SIMP')||strcmp(optMode,'SIMP')||...
         strcmp(OldoptMode,'SEXT')||strcmp(optMode,'SEXT')) && ...
         not(strcmp(OldoptMode,optMode)) )
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
        if(strcmp(contf,'Yes')&&isfield(MOGAResults,'DVLins'))
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

% sets fixed linear lattice parameters in case of pure sextupole
% optimization case 
%
if (strcmp(optMode,'SEXT'))
  LatticeOptData.HACHRO = setLins(1,LatticeOptData.HACHRO,LatticeOptData,DVLins);
  LatticeOptData.ACHRO  = setLins(2,LatticeOptData.ACHRO,LatticeOptData,DVLins);
  LatticeOptData.UC     = setLins(3,LatticeOptData.UC,LatticeOptData,DVLins);
  LatticeOptData.RING   = setLins(5,LatticeOptData.RING,LatticeOptData,DVLins);
end

Constfnct=@(x)stabcon(x,LatticeOptData);

tic;
[x,fval,exitflag,output,population,score] = gamultiobj...
(LattOptfnct,nvars,[],[],[],[],lb,ub,Constfnct,options);

%Collects results in output structure
MOGAResults.LatticeOptData=LatticeOptData;
MOGAResults.options=options;
MOGAResults.opffnct=optfnct;
MOGAResults.fitvars=fitvars;
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
if(strcmp(optMode,'SEXT'))
    MOGAResults.DVLins=DVLins;
end

toc;
nfitvars = size(fitvars,1);
%
% complements lattice data for Pareto Front Indviduals
%
if (strcmp(compf,'Y'))
   fprintf('Colllecting additional ParetoFront data \n'); 
   pfsize=size(MOGAResults.ParetoFront,1); 
   tic;
   for i=1:pfsize
       try
        rpar=ExMOGA(MOGAResults,-i,'N',1,'N', 'N', 'N', 'N', 'N', [], [],[],'N','N',LatticeOptData.DAoptions);
        MOGAResults.ParetoFront(i,nvars+nfitvars+1)=rpar.tunesuc0(1);
        MOGAResults.ParetoFront(i,nvars+nfitvars+2)=rpar.tunesuc0(2);
        MOGAResults.ParetoFront(i,nvars+nfitvars+3)=rpar.beta0(1);
        MOGAResults.ParetoFront(i,nvars+nfitvars+4)=rpar.beta0(2);
        MOGAResults.ParetoFront(i,nvars+nfitvars+5)=rpar.etax;
        MOGAResults.ParetoFront(i,nvars+nfitvars+6)=rpar.SD;
        MOGAResults.ParetoFront(i,nvars+nfitvars+7)=rpar.SFI;
        MOGAResults.ParetoFront(i,nvars+nfitvars+8)=i;
       catch ME
          fprintf('Error getting linear optics for individual n. %4d \n', i);
          fprintf('Error message was:%s \n',ME.message);
       end
   end
   toc;
end
if(strcmp(savef,'Y'))
    fprintf('Saving file %s \n', filename);
    filename=strcat('MOGA_Scans/',filename);
    save(filename,'MOGAResults');
end

fprintf('----- \n');
fprintf('MOGA Scan ended on %s  \n', datetime);
fprintf('----- \n');

end


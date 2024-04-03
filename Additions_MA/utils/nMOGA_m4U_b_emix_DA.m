% --------------------------------------------------------
% MOGA optimisation of emittance and dynamic aperture
% ---------------------------------------------------------
function [X,fval,exitflag,output,MOGAResults] = nMOGA_m4U_b_emix_DA(X0)%, besti)
if nargin<1
    X0 = [];
end

%MOGAdir = '/home/magsjo/ControlRoom/at/atmat/pubtools/MOGAissues';
%persistent ri; 

%ri = m4U_240114_b01_02_03_02__grd_segmented;
ri = m4U_240314_b01_02_04_03;
if ~isempty(X0)
    ri  = alter_m4U_lattice(X0, ri, 'B'); % if new intial parameters alter the lattice
else
    X0 = retrieve_initial_parameters(ri); % if not retrieve parameters and keep the lattice
end

% % %  dx     = [0.09,     0.09,    0.09,     0.09,     0.09,     0.09,     0.09,    0.6*4,  0.6*4,  0.6*4,  0.6*4,  0.6*4];
% % % dx     = [0.10, 0.10, 0.10, 0.10, 0.10, 0.10, ... % QUAD
% % %           0.10, 0.10, ...                         % RB quad   
% % %           2.00, 2.00, ...                         % BB quad 
% % %           2.00, 2.00, 2.00, 2.00, 2.00, ...       % SEXT   
% % %           2.00, 2.00, 2.00]*1e-3;                % OCT  
dx     = [2.00, 2.00, 2.00, 2.00, 2.00, 2.00, ... % QUAD
          2.00, 2.00, ...                         % RB quad   
          4.00, 4.00, ...                         % BB quad 
          8.00, 8.00, 8.00, 8.00, 8.00, ...       % SEXT   
          100.00, 100.00, 100.00]*2.5e-3;         % OCT  
nvars  = length(dx);
sizPop   = 200; % 400;  % 5000
MaxGen   = 5; %15; % 20
FuncTol  = 1e-3;
ConstTol = 1e-1;

if (size(X0,2)==nvars)
    lb = X0-abs(X0).*dx/1;
    ub = X0+abs(X0).*dx/1;
    %fprintf('lower bounds reset to \r\n');
    %fprintf('upper bounds reset to \r\n');
end

if 1==0 %old way
    iniPop = 2*(rand(sizPop,nvars)-0.5).*dx+X0;
else % a la Pedro
    %cd(MOGAdir);
    [MOGAfilename, MOGAdir]=uigetfile('../MOGA_*.mat');
    if (ischar(MOGAfilename))
        load([MOGAdir, MOGAfilename],'MOGAResults');
    else
        fprintf('No file chosen. Aborting... \n');
        return
    end
    OldPopulation     = MOGAResults.population;
    OldX              = MOGAResults.X;
    OldPopSize        = MOGAResults.sizPop;
    Oldnvars          = MOGAResults.nvars;
    mvars             = min(nvars,Oldnvars);
    mPopSize          = min(OldPopSize,sizPop);
    iniPop = zeros(sizPop,nvars);
    for i=1:sizPop
          for j=1:nvars
               iniPop(i,j) = (lb(j)+ub(j))/2;
          end
    end
    iniPop(1:mPopSize,1:mvars)=OldPopulation(1:mPopSize,1:mvars);
    %iniPop(1:mPopSize,1:mvars)=OldX(1:mPopSize,1:mvars);
   
end

% iniPop = [iniPop; X0]; 

%sizIni = 100;
%iniPop = 2*(rand(sizIni,12)-0.5).*delta/100+X0;
%iniPop =  repmat(X0,sizPop,1);


% iniPop = 2*(rand(5000-size(X0,1),12)-0.5).*delta+X0(besti,:);
% iniPop = [iniPop; X0]; 

%iniPop = X0; 
controlgcp = 0; 
usepar = false;
if controlgcp==1
   c = parcluster;
   % next line is optional at MAX IV
   c.AdditionalProperties.AccountName = 'any-virtual-account-name';
   % 6 hour walltime
   c.AdditionalProperties.WallTime = '70:00:00';
   % hyperthreading enabled
   c.saveProfile
   c.NumThreads = 1; %2

   poolobj = gcp('nocreate'); % get the pool object, and do it avoiding creating a new one.
   if isempty(poolobj) % check if there is not a pool.
    %poolsize = 1;
   else
    %poolsize = poolobj.NumWorkers;
    delete( gcp('nocreate')); % delete the current pool object.
   end
   parpool('aurora R2022a',96); % 128 %48 %96

   % % % %    c.NumThreads = 1;
% % % %    parpool('local',12) %48
   pp = gcp; 
   %% usepar = true;
end



%
% start with the initial options
%
% options = optimoptions('gamultiobj');
options = optimoptions('gamultiobj','InitialPopulationRange',[lb;ub]);

%% Modify options setting
% options = optimoptions(options,'CreationFcn', @gacreationnonlinearfeasible);
% options = optimoptions(options,'MutationFcn', @mutationadaptfeasible);
% options = optimoptions(options,'CrossoverFcn',@crossoverintermediate);

% options = optimoptions(options,'PlotFcn',@gaplotpareto);
options = optimoptions(options,'PlotInterval',1);
%options = optimoptions(options,'display', 'iter');
options = optimoptions(options,'PopulationSize',sizPop);
options = optimoptions(options,'InitialPopulation',iniPop);
%options = optimoptions(options,'InitialPopulationMatrix',iniPop);
% options = optimoptions(options,'InitialPop',iniPop);
options = optimoptions(options,'UseParallel',usepar);
options = optimoptions(options,'UseVectorized', false);
options = optimoptions(options,'OutputFcns', @MOGAoutputfcn);
options = optimoptions(options,'FunctionTolerance',FuncTol);%1e-3);
options = optimoptions(options,'MaxGenerations',MaxGen);%40
options = optimoptions(options,'ConstraintTolerance',100);   %ConstTol);%1e-1);
options = optimoptions(options,'PlotFcn',{@gaplotpareto,@gaplotrankhist});


% myobj      = @(X)OBJfcn_m4U(X, ri); % object function
% mycon      = @(X)CONfcn_m4U(X, ri); % constraint function
myobj      = @(x) nOBJfcn_m4U(x,ri); % object function
mycon      = @(x) nCONfcn_m4U(x,ri); % constraint function
nvars      = length(X0);
%[X,fval,exitflag,output,population,score] = gamultiobj(myobj,nvars,[],[],[],[],[], [], mycon, options);
[X,fval,exitflag,output,population,score] = gamultiobj(myobj,nvars,[],[],[],[],lb, ub, mycon, options);
%[X,fval,exitflag,output,population,score] = gamultiobj(myobj,nvars,[],[],[],[],lb, ub, options);

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


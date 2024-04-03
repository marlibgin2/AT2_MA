function [state, options,optchanged] = MOGAoutputfcn(options,state,flag)
%displays the function eval value at each iteration. You can change this
disp(state.FunEval);
optchanged = false;
switch flag
 case 'init'
        disp('Starting the algorithm');
    case {'iter','interrupt'}
        disp('Iterating ...')
        adesso = datestr(now,'ddmmyyyy_HH-MM-SS');
        save(['nMOGAgen_' adesso '.mat'],'state','options')
    case 'done'
        disp('Performing final task');
        adesso = datestr(now,'ddmmyyyy_HH-MM-SS');
        save(['nMOGAfin_' adesso '.mat'],'state','options')
end
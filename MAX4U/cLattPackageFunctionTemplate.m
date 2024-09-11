function OutStruct = cLattPackageFunctionTemplate(varargin)
% This is a template for new functions to add functionality to the cLatt package
% 
% Explanation of what the function does
%% Inputs
% Mandatory arguments
%
% Optional  arguments
% verbose              :defines level of verbose output, default=0, i.e. no output
% Optional flags
% 
%% Outputs
% OutStruct.inputs : echoes inputs
%
% OutStruct.outputs.out1 : description
% OutStruct.outputs.out1 : decription
%
%% Usage Examples
% 

% note: the blank line above is important !
%% History
% NN 2024/xx/xx : first version
% MM date       : updates
%
%% Preamble
% any required initializations


%% Input argument parsing
LattSt               = getargs(varargin,[]); % Example reading of maadatory (positional) argument

lattname             = getoption(varargin,'lattname',''); % Example reading of optional argument

verboselevel         = getoption(varargin,'verbose',0);

basicf               = any(strcmpi(varargin,'basic')); % Example optional flag

%% Main function body, section 1
tstart=tic;

%% section 2

calcSomething(1,2,'verbose',verboselevel-1) % when calling other functions pass the verboselevel argument

telapsed=toc(tstart);
%% Collects output structure

OutStruct.inputs.field1 = 1; % echos the nmnlut
OutStruct.inputs.field1 = 2;

OutStruct.outputs.field1 = 0;
OutStruct.telapsed = telapsed;




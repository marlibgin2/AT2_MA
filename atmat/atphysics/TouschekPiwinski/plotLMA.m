function plotLMA(varargin)
% plots the Local Momenttum Aperture
%   
%% Inputs
% Mandatory argument
% LMA : structure generated by calcLMA
% 
% Optional arguments
% verbose : defines level of verbose output, default=0, i.e. no output
%
%% Usage examples
% plotLMA(LMA);


%% History
% PFT 2024/06/15, first version
%
%% Input argument parsing
%
LMA           = getargs(varargin,[]);
verboselevel  = getoption(varargin,'verbose',0);


%% Plots LMA
Spos  = LMA.outputs.Spos;
map_l = LMA.outputs.map_l;
map_h = LMA.outputs.map_h;

figure;plot(Spos, map_l*100, '-o');hold;plot(Spos,map_h*100,'o-');
xlabel('S[m]');
ylabel('Local Momentum Aperture [%]');
grid on;
ylim([-15,15]);




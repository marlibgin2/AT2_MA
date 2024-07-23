function plotLMA(varargin)
% plots the Local Momenttum Aperture
%   
%% Inputs
% Mandatory argument
% LMA : structure generated by calcLMA
% 
% Optional arguments
% dpmaxplot : maximum of energy deviation axis, default =  LMA.deltalimit
% dpminplot : minimum of energy deviation axis, default = -LMA.deltalimit
% plottitle : title string
% verbose : defines level of verbose output, default=0, i.e. no output
%
%% Usage examples
% plotLMA(LMA);
% plotLMA(LMA,'dpmaxplot', 0.30, 'dpminplot', -0.3);

%% History
% PFT 2024/06/15, first version
% PFT 2024/07/03 added vertical scale control
% PFT 2024/07/15 added optional plot title
%% Input argument parsing
%
LMA           = getargs(varargin,[]);
verboselevel  = getoption(varargin,'verbose',0);
dpmaxplot     = getoption(varargin,'dpmaxplot', LMA.outputs.MAoptions.deltalimit);
dpminplot     = getoption(varargin,'dpminplot',-LMA.outputs.MAoptions.deltalimit);
plottitle     = getoption(varargin,'plottitle','');

%% Plots LMA
Spos  = LMA.outputs.Spos;
map_l = LMA.outputs.map_l;
map_h = LMA.outputs.map_h;

figure;plot(Spos, map_l*100, '-o');hold on;plot(Spos,map_h*100,'o-');
xlabel('S[m]');
ylabel('Local Momentum Aperture [%]');
grid on;
ylim([dpminplot*100,dpmaxplot*100]);
title(strcat(plottitle, {' LMA without errors'}));




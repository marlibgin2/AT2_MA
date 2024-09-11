function plotOrbit(sBPM, orb, varargin)
% Plots the closed orbit
%% Inputs
% Mandatory arguments
% sBPM: (1XN) array of BPM positions [m]
% orb : (6XN) array of Closed Orbit coordinates ([m] [rad] [m] [rad] [] [m])
%
% Optional arguments
%
% S0   : min X axis plot range, default = 0.0
% Smax : max X Axis plot range, defaul = max(SBPm)
%
%% Input arumenta parsing
S0   = getoption(varargin,'S0',0.0);
Smax = getoption(varargin,'Smax',max(sBPM));
%% Plots the orbit
figure; plot(sBPM,1e3*orb([1 3],:)); xlim([S0 Smax]);
            xlabel('S[m]'); ylabel('x,y [mm]');grid;legend('X','Y');
 
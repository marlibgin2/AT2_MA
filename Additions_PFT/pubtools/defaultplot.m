function plotdata=defaultplot(lindata,ring,dpp,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 beta=cat(1,lindata.beta);                     % left axis
 plotdata(1).values=beta;
 plotdata(1).labels={'\beta_x','\beta_z'};
 plotdata(1).axislabel='\beta [m]';
 dispersion=cat(2,lindata.Dispersion)';        % right axis
 plotdata(2).values=dispersion(:,1);
 plotdata(2).labels={'\eta_x'};
 plotdata(2).axislabel='dispersion [m]';
end


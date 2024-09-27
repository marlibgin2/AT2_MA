function phandles=plotTLTdist(varargin)
% plots distribution of Touschek LIfetimes 
%% Inputs
% Mandatory argument
% TLdist : structure generated by calcTLTdist
% 
% optional arguments
% dpmaxplot : maximum of energy deviation axis, default =  LMA.deltalimit
% dpminplot : minimum of energy deviation axis, default = -LMA.deltalimit
% plottitle : title string
%
% Optional flags
% plotorbrms: plots the rms of the uncorrected and corrected orbits for all
%             seeds

%
%% Usage examples
% plotTLTdist(TLdist);
% plotTLTdist(TLdist,'plotorbrms');
% plotLMAdist(TLdist,'plotorbrms','verbose',1);
%
% see also calcTLTdist

%% History
% PFT 2024/07/05, first version
% PFT 2024/07/15: added optional title string

%% Input argument parsing
%
TLdist       = getargs(varargin,[]);
plotorbrmsf  = any(strcmpi(varargin,'plotorbrms'));
plottitle    = getoption(varargin,'plottitle','');
verboselevel = getoption(varargin,'verbose',0);

%% Plots Touschek lifetime 
nhandles=0;
phandles={};

nseeds=TLdist.inputs.nseeds;
figure; plot(0:nseeds,TLdist.outputs.TLs/3600, '-o')
% -----------------------------------------
% MA 27092024
% calculate and print average TLT +/- error
% -----------------------------------------
aTLT = mean(TLdist.outputs.TLs(2:end)/3600);
sTLT = std(TLdist.outputs.TLs(2:end)/3600); 
ypl   = get(gca,'YLim');
xpl   = get(gca,'XLim');
text(mean(xpl)-2,mean(ypl),['<TLT> = ' num2str(aTLT,3) ' +/- ' num2str(sTLT,2) ' (hr)'])
xlabel('Seed');ylabel('Touschek LIfetime [hr]');
grid on;
title(plottitle);
nhandles=nhandles+1;
phandles{1}=gcf;

if (plotorbrmsf)
   figure; 
   plot(TLdist.outputs.LMAdist.outputs.orb0_stds(1,2:end)*1000,'-o');
   xlabel('seed #');
   ylabel('x/y[mm]');
   hold on;
   plot(TLdist.outputs.LMAdist.outputs.orb0_stds(3,2:end)*1000,'-o');
   legend({'X','Y'});
   title(strcat(plottiel,'rms orbit before correction'));
   nhandles=nhandles+1;
   phandles{nhandles}=gcf;


   figure; 
   plot(TLdist.outputs.LMAdist.outputs.orb_stds(1,2:end)*1000,'-o');
   xlabel('seed #');
   ylabel('x/y[mm]');
   hold;
   plot(TLdist.outputs.LMAdist.outputs.orb_stds(3,2:end)*1000,'-o');
   legend({'X','Y'});
   title(strcat(plottitle,'rms orbit after correction'));
   nhandles=nhandles+1;
   phandles{nhandles}=gcf;
end



 



function phandles=plotTuneMap(tunemap,varargin)
% Plots tune maps 
%% Inputs  
% Mandatory arguments
% tunemap: structure generated by the calcTuneMap function
%
% Optional arguments
%
% note: defaults are taken from tunemap input structure
%
% plotmode: 'abs' : plots full tune value (including integer part)
%           'rel' : plots tune variations with respect to small amplitude
%                   tunes.
%            note: plotmode is onlöyrelevane for plottypes :x,y,xy
% plottype: 
%   'x'       : ADTS along horizontal axis (default)
%   'y'       : ADTS along vertical axis
%   'xy'      : both ADTS
%   'td'      : same as 'xy', but on a tune diagram
%   'gridtd'  : grid of points on a tune diagram
%   'gridxy'  : grid of points as a series of lines vs 
%               x for a chosen set of y's
%   'gridyx'  : grid of points as a series of lines vs 
%               y for a chosen set of x's
%   'gridxdp' : grid of points as a series of lines vs 
%               x for a chosen set of dp's
%   'gridydp' : grid of points as a series of lines vs 
%               y for a chosen set of dp's
%   'griddpx' : grid of points as a series of lines vs 
%               dp for a chosen set of x's
%   'griddpy' : grid of points as a series of lines vs 
%               dp for a chosen set of y's
%   'difxy' : tune diffusion plot on xy plane
%   'difxdp': tune diffusion plot on xdp plane
%   'difydp': tune diffusion plot on ydp plane
%   'fmxy'  : xy tune diffusion plot on tune diagram (frequency map)
%   'fmxdp' : xdp tune diffusion plot on tune diagram (frequency map) 
%   'fmydp' : xdp tune diffusion plot on tune diagram (frequency map)
%   'chro'  : tunes vs energy deviation
%   'chrotd': chromatic tune footprint on a tune diagram
%
% resorder: resonance order for tune diagram, default = 3
% qxrange=[qxmin qymin]: horizontal plot range in tune diagram,  default =[0 1]
% qyrange=[qymin qymax]: vertical plot range in tune diagram, default= [0 1]
% xminplot: minimum horizontal coordinate in tune plots, default from tunemap structure;
% xmaxplot: maximum horizontal coordinate in tune plots, default from tunemap structure;
% yminplot: minimum vertical coordinate in tune plots, default from tunemap structure;
% ymaxplot: maximum vertical coordinate in tune plots, default from tunemap structure;
% dpminplot: minimum momentum deviation in tune plots, deafult = dpmin
% dpmaxplot: maximum momentum deviation in tune plots, default = dpmax
% xminplot_dm: minimum horizontal coordinate in diffusion maps, default from tunemap structure;
% xmaxplot_dm: maximum horizontal coordinate in diffusion maps, default from tunemap structure;
% yminplot_dm: minimum vertical coordinate in diffusion maps, default from tunemap structure;
% ymaxplot_dm: maximum vertical coordinate in diffusion maps, default from tunemap structure;
% dpminplot_dm: minimum momentum deviation in diffusion maps, deafult = dpmin
% dpmaxplot_dm: maximum momentum deviation in diffusion maps, default = dpmax
% x0  : (1xN) array of horizontal values for line plots vs y or dp
% y0  : (1xN) array of horizontal values for line plots vs x or dp
% dp0 : (1XN) array of energy deviation values  for line plots vs x or y
%
% caxrange= [cmin cmax]: color axis range, default from input structure
% dqx     : horizontal half width of square in tune space for FMA plots,
%           default = 0.001
% dqy     : vertical half width of square in tune space for FMA plots,
%           default = 0.001 
%
% plottitle   : string to be added to the plot title, default = '' 
%
% Optional flags
% rate        : selects diffusion rate instead of diffusion
%
%% Outputs
% phandles = cell array with handles to created plots
%
%% Usage examples
% plotTuneMap(tunemap,'plottype','x','plotmode','rel');
% plotTuneMap(tunemap,'plottype','td','qxrange',[0.2 0.3],'resorder',3);
% plotTuneMap(tunemap,'plottype','fmxy','qxrange',[0.2 0.3],'resorder',3);
% plotTuneMap(tunemap,'plottype','gridxy','y0',0.0);

%% History
% PFT 2024/04/27: first version, based on plotADTS
% PFT 2024/05/03: added chromatic tunemap plots
% PFT 2024/05/05: added tune diffusion plot in (x,y) plane
% PFT 2024/05/08: added tune diffusion plots in (x,dp) and (y,dp) planes
% PFT 2024/05/10: added grid plots for tune variations
% PFT 2024/07/05: updated to handle TMoptions structure from calcTuneMap
% PFT 2024/07/06: changed hold statement to avoid clutter
% PFT 2024/07/12: separated arguments for plot ragtes of tune maps and
%                 diffusion maps
% PFT 2024/07/13: added string to plot title
% PFT 2024/07/15: added output handles, reformatted titles
% PFT 2024/07/16: incoporrated new plot_net function that can handle
%                 the integer part of the tune

%% Input argument parsing
plotmode   = getoption(varargin,'plotmode',tunemap.inputs.plotargs.plotmode);
plottype   = getoption(varargin,'plottype',tunemap.inputs.plotargs.plottype);
resorder   = getoption(varargin,'resorder',tunemap.inputs.plotargs.resorder);
qxrange    = getoption(varargin,'qxrange',tunemap.inputs.plotargs.qxrange);
qyrange    = getoption(varargin,'qyrange',tunemap.inputs.plotargs.qyrange);
caxrange   = getoption(varargin,'caxrange',tunemap.inputs.plotargs.caxrange);
plottitle  = getoption(varargin,'plottitle','');
dqx        = getoption(varargin,'dqx',tunemap.inputs.plotargs.dqx);
dqy        = getoption(varargin,'dqy',tunemap.inputs.plotargs.dqy);
x0         = getoption(varargin,'x0',tunemap.inputs.plotargs.x0);
y0         = getoption(varargin,'y0',tunemap.inputs.plotargs.y0);
dp0        = getoption(varargin,'dp0',tunemap.inputs.plotargs.dp0);
ratef      = any(strcmpi(varargin,'rate'));

xmin       = tunemap.outputs.TMoptions.xmin;
xmax       = tunemap.outputs.TMoptions.xmax;
ymin       = tunemap.outputs.TMoptions.ymin;
ymax       = tunemap.outputs.TMoptions.ymax;
dpmin      = tunemap.outputs.TMoptions.dpmin;
dpmax      = tunemap.outputs.TMoptions.dpmax;
xmin_dm    = tunemap.outputs.TMoptions.xmin_dm;
xmax_dm    = tunemap.outputs.TMoptions.xmax_dm;
ymin_dm    = tunemap.outputs.TMoptions.ymin_dm;
ymax_dm    = tunemap.outputs.TMoptions.ymax_dm;
dpmin_dm   = tunemap.outputs.TMoptions.dpmin_dm;
dpmax_dm   = tunemap.outputs.TMoptions.dpmax_dm;

npx        = tunemap.outputs.TMoptions.npx;
npy        = tunemap.outputs.TMoptions.npy;
npd        = tunemap.outputs.TMoptions.npd;
amplx      = tunemap.outputs.amplx;
amply      = tunemap.outputs.amply;
axgridxy   = tunemap.outputs.axgridxy;
aygridxy   = tunemap.outputs.aygridxy;
axgridxdp  = tunemap.outputs.axgridxdp;
aygridydp  = tunemap.outputs.aygridydp;
dpgridxdp  = tunemap.outputs.dpgridxdp;
dpgridydp  = tunemap.outputs.dpgridydp;
dps        = tunemap.outputs.dps;

xminplot   = getoption(varargin,'xminplot',xmin);
xmaxplot   = getoption(varargin,'xmaxplot',xmax);
yminplot   = getoption(varargin,'yminplot',ymin);
ymaxplot   = getoption(varargin,'ymaxplot',ymax);
dpminplot  = getoption(varargin,'dpminplot',dpmin);
dpmaxplot  = getoption(varargin,'dpmaxplot',dpmax);

xminplot_dm   = getoption(varargin,'xminplot_dm',xmin_dm);
xmaxplot_dm   = getoption(varargin,'xmaxplot_dm',xmax_dm);
yminplot_dm   = getoption(varargin,'yminplot_dm',ymin_dm);
ymaxplot_dm   = getoption(varargin,'ymaxplot_dm',ymax_dm);
dpminplot_dm  = getoption(varargin,'dpminplot_dm',dpmin_dm);
dpmaxplot_dm  = getoption(varargin,'dpmaxplot_dm',dpmax_dm);


%% Plots Tune Map
nhandles=0;
phandles={};
switch plotmode
    case 'abs'
       Qxxplot   = tunemap.outputs.Qxx;
       Qyxplot   = tunemap.outputs.Qyx;
       Qxyplot   = tunemap.outputs.Qxy;
       Qyyplot   = tunemap.outputs.Qyy;
       Qxdpplot  = tunemap.outputs.Qxdp;
       Qydpplot  = tunemap.outputs.Qydp;
       Qxgridxyplot = tunemap.outputs.Qxgridxy;
       Qygridxyplot = tunemap.outputs.Qygridxy;
       Qxgridxdpplot = tunemap.outputs.Qxgridxdp;
       Qygridxdpplot = tunemap.outputs.Qygridxdp;
       Qxgridydpplot = tunemap.outputs.Qxgridydp;
       Qygridydpplot = tunemap.outputs.Qygridydp;

    case 'rel'
        Qxxplot  = tunemap.outputs.dQxx;
        Qyxplot  = tunemap.outputs.dQyx;
        Qxyplot  = tunemap.outputs.dQxy;
        Qyyplot  = tunemap.outputs.dQyy;
        Qxdpplot = tunemap.outputs.dQxdp;
        Qydpplot = tunemap.outputs.dQydp;
        Qxgridxyplot = tunemap.outputs.dQxgridxy;
        Qygridxyplot = tunemap.outputs.dQygridxy;
        Qxgridxdpplot = tunemap.outputs.dQxgridxdp;
        Qygridxdpplot = tunemap.outputs.dQygridxdp;
        Qxgridydpplot = tunemap.outputs.dQxgridydp;
        Qygridydpplot = tunemap.outputs.dQygridydp;

    otherwise
        fprintf('%s Error in plotTuneMap. Unknown plot mode : %s \n', ...
                   datetime, plotmode);
        return
end

switch plottype

    case {'x';'X'}
        if not(isempty(amplx)||isempty(Qxxplot)||isempty(Qyxplot))
            figure; xlim([xminplot*1000,xmaxplot*1000]);
            plot(amplx*1000,Qxxplot,'-ok');xlabel('X[mm]');
            hold on;
            if (strcmpi(plotmode,'abs'))
                ylabel('Qx');
                yyaxis right; 
            else
                ylabel('dQ');
            end
            plot(amplx*1000,Qyxplot,'-or');
            if (strcmpi(plotmode,'abs'))
                ylabel('Qy');
                legend({'Qx';'Qy'});
            else
                legend({'dQx';'dQy'});
            end
            grid on;
            title(plottitle);
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qxx or Qyx data not available in tunemap structure for plottype %s\n', datetime, plottype);
        end

    case {'y';'Y'}
        if not(isempty(amply)||isempty(Qxyplot)||isempty(Qyyplot))
            figure; xlim([yminplot,ymaxplot]*1000); 
            plot(amply*1000,Qxyplot,'-ok');xlabel('Y[mm]');
            hold on;
            if (strcmpi(plotmode,'abs'))
                ylabel('Qx');
                yyaxis right; 
            else
                ylabel('dQ');
            end
            plot(amply*1000,Qyyplot,'-or');
            if (strcmpi(plotmode,'abs'))
                ylabel('Qy');
                legend({'Qx';'Qy'});
            else
                legend({'dQx';'dQy'});
            end
            grid on;
            title(plottitle);
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qxy or Qyy data not available in tunemap structure for plottype %s\n', datetime, plottype);
        end

    case {'xy';'XY';'xY';'Xy'}
        if not(isempty(amplx)||isempty(Qxxplot)||isempty(Qyxplot))
            figure; xlim([xminplot,xmaxplot]*1000);
            plot(amplx*1000,Qxxplot,'-ok');xlabel('X[mm]');
            hold on;
            if (strcmpi(plotmode,'abs'))
                ylabel('Qx');
                yyaxis right; 
            else
                ylabel('dQ');
            end
            plot(amplx*1000,Qyxplot,'-or');
            if (strcmpi(plotmode,'abs'))
                ylabel('Qy');
                legend({'Qx';'Qy'});
            else
                legend({'dQx';'dQy'});
            end
            grid on;
            title(plottitle);
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
              fprintf('%s Error in plotTuneMap: Qxx or Qyx data not available in tunemap structure for plottype %s\n', datetime, plottype);
        end

        if not(isempty(amply)||isempty(Qxyplot)||isempty(Qyyplot))
            figure; xlim([yminplot,ymaxplot]*1000); 
            plot(amply*1000,Qxyplot,'-ok');xlabel('Y[mm]');
            hold on;
            if (strcmpi(plotmode,'abs'))
                ylabel('Qx');
                yyaxis right; 
            else
                ylabel('dQ');
            end
            plot(amply*1000,Qyyplot,'-or');
            if (strcmpi(plotmode,'abs'))
                ylabel('Qy');
                legend({'Qx';'Qy'});
            else
                legend({'dQx';'dQy'});
            end
            grid on;
            title(plottitle);
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
              fprintf('%s Error in plotTuneMap: Qxy or Qyy data not available in tunemap structure for plottype %s\n', datetime, plottype);
        end

    case {'td';'TD'}
        Qxx=tunemap.outputs.Qxx;
        Qyx=tunemap.outputs.Qyx;
        Qxy=tunemap.outputs.Qxy;
        Qyy=tunemap.outputs.Qyy;

        if not(isempty(Qxx)||isempty(Qyx)||isempty(Qyx)||isempty(Qyy))
            figure;plot(Qxx,Qyx,'ok');hold on;
            plot(Qxy,Qyy,'or');
            plot_net(resorder,qxrange(1),qxrange(2),...
                          qyrange(1),qyrange(2));
                          
            xlabel('qx');ylabel('qy');
            title(strcat(plottitle,{ 'Res order = '},num2str(resorder)));
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qxx, Qxy, Qyx or Qyy not available for plottype %s \n', datetime, plottype);
        end

    case {'gridtd';'GRIDTD'}
        Qxgridxy=tunemap.outputs.Qxgridxy;
        Qygridxy=tunemap.outputs.Qygridxy;

        if (not(isempty(Qxgridxy)||isempty(Qygridxy)))
            figure;plot(Qxgridxy,Qygridxy,'ok');hold on;
            plot_net(resorder,qxrange(1),qxrange(2),...
                          qyrange(1),qyrange(2));
                          
            xlabel('Qx');ylabel('Qy');
            title(strcat(plottitle,{' Res order = '},num2str(resorder)));
            grid on;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qxgrid or Qygrid not available for plottype %s \n', datetime, plottype);
        end

    case {'gridxy';'GRIDXY'}
        if (not(isempty(Qxgridxyplot)||isempty(Qygridxyplot)||isempty(amplx)||isempty(amply)))
            if (ischar(y0))
                plotslice(amplx*1000,amply*1000,Qxgridxyplot,'x',y0);
            else
                plotslice(amplx*1000,amply*1000,Qxgridxyplot,'x',y0*1000);
            end
            xlim([xminplot xmaxplot]*1000);
            xlabel('X[mm]');ylabel('qx');
            title(strcat(plottitle,{' Qx vs x at fixed y'}));
            grid on;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
            if (ischar(y0))
                plotslice(amplx*1000,amply*1000,Qygridxyplot,'x',y0);
            else
                plotslice(amplx*1000,amply*1000,Qygridxyplot,'x',y0*1000);
            end
            xlim([xminplot xmaxplot]*1000);
            xlabel('X[mm]');ylabel('qy');
            title(strcat(plottitle,{' Qy vs x at fixed y'}));
            grid on;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qxgridxy or Qygridxy not available for plottype %s \n', datetime, plottype);
        end

    case {'gridyx';'GRIDYX'}
        if (not(isempty(Qxgridxyplot)||isempty(Qygridxyplot)||isempty(amplx)||isempty(amply)))
            if (ischar(x0))
                plotslice(amplx*1000,amply*1000,Qxgridxyplot,'y',x0);
            else
                plotslice(amplx*1000,amply*1000,Qxgridxyplot,'y',x0*1000);
            end
            xlim([yminplot ymaxplot]*1000);
            xlabel('Y[mm]');ylabel('qx');
            title(strcat(plottitle,{' Qx vs y at fixed x'}));
            grid on;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
            if (ischar(x0))
                plotslice(amplx*1000,amply*1000,Qygridxyplot,'y',x0);
            else
                plotslice(amplx*1000,amply*1000,Qygridxyplot,'y',x0*1000);
            end
            xlim([yminplot ymaxplot]*1000);
            xlabel('Y[mm]');ylabel('qy');
            title(strcat(plottile,{' Qy vs y at fixed x'}));
            grid on;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qxgridxy or Qygridxy not available for plottype %s \n', datetime, plottype);
        end

     case {'gridxdp';'GRIDXDP'}
        if (not(isempty(Qxgridxdpplot)||isempty(Qygridxdpplot)||isempty(amplx)||isempty(dps)))
            if (ischar(dp0))
                plotslice(dps*100,amplx*1000,Qxgridxdpplot,'y',dp0);
            else
                plotslice(dps*100,amplx*1000,Qxgridxdpplot,'y',dp0*100);
            end
            xlim([xminplot xmaxplot]*1000);
            xlabel('X[mm]');ylabel('qx');
            title('Qx vs x at fixed dp');
            grid on;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
            if (ischar(dp0))
                plotslice(dps*100,amplx*1000,Qygridxdpplot,'y',dp0);
            else
                plotslice(dps*100,amplx*1000,Qygridxdpplot,'y',dp0*100);
            end
            xlim([xminplot xmaxplot]*1000);
            xlabel('X[mm]');ylabel('qy');
            title(strcat(plottitle,{' Qy vs x at fixed dp'}));
            grid on;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qxgridxdp or Qygridxdp not available for plottype %s \n', datetime, plottype);
        end

     case {'griddpx';'GRIDDPX'}
        if (not(isempty(Qxgridxdpplot)||isempty(Qygridxdpplot)||isempty(amplx)||isempty(dps)))
            if (ischar(x0))
                plotslice(dps*100,amplx*1000,Qxgridxdpplot,'x',x0);
            else
                plotslice(dps*100,amplx*1000,Qxgridxdpplot,'x',x0*1000);
            end
            xlim([dpminplot dpmaxplot]*100);
            xlabel('dp[%]');ylabel('qx');
            title(strcat(plottile,{' Qx vs dp at fixed x}'}));
            grid on;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
            if (ischar(x0))
                plotslice(dps*100,amplx*1000,Qygridxdpplot,'x',x0);
            else
                plotslice(dps*100,amplx*1000,Qygridxdpplot,'x',x0*1000);
            end
            xlim([dpminplot dpmaxplot]*100);
            xlabel('dp[%]');ylabel('qy');
            title(strcat(plottile,{' Qy vs dp at fixed x'}));
            grid on;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qxgridxdp or Qygridxdp not available for plottype %s \n', datetime, plottype);
        end

    case {'gridydp';'GRIDtDP'}
        if (not(isempty(Qxgridydpplot)||isempty(Qygridydpplot)||isempty(amply)||isempty(dps)))
            if (ischar(dp0))
                plotslice(dps*100,amply*1000,Qxgridydpplot,'y',dp0);
            else
                plotslice(dps*100,amply*1000,Qxgridydpplot,'y',dp0*100);
            end

            xlim([yminplot ymaxplot]*1000);
            xlabel('Y[mm]');ylabel('qx');
            title(strcat(plottitle,{' Qx vs y at fixed dp'}));
            grid on;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
            if (ischar(dp0))
                plotslice(dps*100,amply*1000,Qygridydpplot,'y',dp0);
            else
                plotslice(dps*100,amply*1000,Qygridydpplot,'y',dp0*100);
            end
            xlim([yminplot ymaxplot]*1000);
            xlabel('Y[mm]');ylabel('qy');
            title(strcat(plottitle,{' Qy vs y at fixed dp'}));
            grid on;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qxgridydp or Qygridydp not available for plottype %s \n', datetime, plottype);
        end

    case {'griddpy';'GRIDDPY'}
        if (not(isempty(Qxgridydpplot)||isempty(Qygridydpplot)||isempty(amply)||isempty(dps)))
            if (ischar(y0))
                plotslice(dps*100,amply*1000,Qxgridydpplot,'x',y0);
            else
                plotslice(dps*100,amply*1000,Qxgridydpplot,'x',y0*1000);
            end
            xlim([dpminplot dpmaxplot]*100);
            xlabel('dp[%]');ylabel('qx');
            title(strcat(plottile,{' Qx vs dp at fixed y'}));
            grid on;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
            if (ischar(y0))
                plotslice(dps*100,amply*1000,Qygridydpplot,'x',y0);
            else
                plotslice(dps*100,amply*1000,Qygridydpplot,'x',y0*1000);
            end
            xlim([dpminplot dpmaxplot]*100);
            xlabel('dp[%]');ylabel('qy');
            title(strcat(plottile,{' Qy vs dp at fixed y'}));
            grid on;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qxgridydp or Qygridydp not available for plottype %s \n', datetime, plottype);
        end

    case {'difxy';'DIFXY'}
        if (ratef)
            Qdifxy   = tunemap.outputs.Qdifxyra;
        else
            Qdifxy   = tunemap.outputs.Qdifxy;
        end
        
        if(not(isempty(Qdifxy)))
            Qdifxymat= reshape(Qdifxy,npy,npx);
        
            figure;
            h=imagesc([min(amplx),max(amplx)]*1000,[min(amply),max(amply)]*1000,...
               Qdifxymat);
            set(h,'alphadata',~isnan(Qdifxymat));
            ax=gca; ax.YDir='normal';
            xlim([xminplot_dm xmaxplot_dm]*1000);ylim([yminplot_dm ymaxplot_dm]*1000);
                 xlabel('X[mm]');ylabel('Y[mm]');
            grid on;
            colormap('jet');
            clim(caxrange);
            shading flat;
            if (ratef)
                title(strcat(plottitle,{' Tune Diffusion Rate'}));
            else
                title(strcat(plottitle,{' Tune Diffusion'}));
            end
            colorbar;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qdifxy not available for plottype %s \n', datetime, plottype);
        end


    case {'difxdp';'DIFXDP'}
        if (ratef)
            Qdifxdp   = tunemap.outputs.Qdifxdpra;
        else
            Qdifxdp   = tunemap.outputs.Qdifxdp;
        end
        
        if (not(isempty(Qdifxdp)))
            Qdifxdpmat= reshape(Qdifxdp,npx,npd);
        
            figure;
            h=imagesc([min(dps),max(dps)]*100,[min(amplx),max(amplx)]*1000,...
               Qdifxdpmat);
            set(h,'alphadata',~isnan(Qdifxdpmat));
            ax=gca; ax.YDir='normal';
            xlim([dpminplot_dm dpmaxplot_dm]*100);ylim([xminplot_dm xmaxplot_dm]*1000);
                 xlabel('dp[%]');ylabel('X[mm]');
            grid on;
            colormap('jet');
            clim(caxrange);
            if (ratef)
                title(strcat(plottitle,{' Tune Diffusion Rate'}));
            else
                title(strcat(plottitle,{' Tune Diffusion'}));
            end
            colorbar;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qdifxdp not available for plottype %s \n', datetime, plottype);
        end


    case {'difydp';'DIFYDP'}
        if (ratef)
            Qdifydp   = tunemap.outputs.Qdifydpra;
        else
            Qdifydp   = tunemap.outputs.Qdifydp;
        end

        if (not(isempty(Qdifydp)))
            Qdifydpmat= reshape(Qdifydp,npy,npd);
        
            figure;
            h=imagesc([min(dps),max(dps)]*100,[min(amply),max(amply)]*1000,...
                      Qdifydpmat);
            set(h,'alphadata',~isnan(Qdifydpmat));
            ax=gca; ax.YDir='normal';
            xlim([dpminplot_dm dpmaxplot_dm]*100);ylim([yminplot_dm ymaxplot_dm]*1000);
                 xlabel('dp[%]');ylabel('Y[mm]');
            grid on;
            colormap('jet');
            clim(caxrange);
            if (ratef)
                title(strcat(plottitle,{' Tune Diffusion Rate'}));
            else
                title(strcat(plottile,{' Tune Diffusion'}));
            end
            colorbar;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qdifydp not available for plottype %s \n', datetime, plottype);
        end


    case {'fmxy';'FMXY'}
        if (ratef)
            Qdifxy   = tunemap.outputs.Qdifxyra;
        else
            Qdifxy   = tunemap.outputs.Qdifxy;
        end

        if(not(isempty(Qdifxy)))
            Qxgridxy=tunemap.outputs.Qxgridxy;
            Qygridxy=tunemap.outputs.Qygridxy;
   
            Qxfmxy(:,1:npx*npy) = [Qxgridxy'-dqx; Qxgridxy'-dqx; ...
                                   Qxgridxy'+dqx; Qxgridxy'+dqx];
            Qyfmxy(:,1:npx*npy) = [Qygridxy'-dqy; Qygridxy'+dqy; ...
                                   Qygridxy'+dqy; Qygridxy'-dqy];
            figure;
            fill(Qxfmxy,Qyfmxy,Qdifxy);hold on;
            xlabel('qx');ylabel('qy');
            plot_net(resorder,qxrange(1),qxrange(2),...
                     qyrange(1),qyrange(2));
            grid on;
            colormap('jet');
            shading flat;
            clim(caxrange);
            if (ratef)
                title(strcat(plottitle,{' Tune Diffusion Rate:(X,Y),Res order='},num2str(resorder)));
            else
                title(strcat(plottitle, {' Tune Diffusion:(X,Y), Res order ='},num2str(resorder)));
            end
            colorbar;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qdifxy not available for plottype %s \n', datetime, plottype);
        end

    case {'fmxdp';'FMXDP'}
        if (ratef)
            Qdifxdp   = tunemap.outputs.Qdifxdpra;
        else
            Qdifxdp   = tunemap.outputs.Qdifxdp;
        end

        if (not(isempty(Qdifxdp)))    
            Qxgridxdp=tunemap.outputs.Qxgridxdp;
            Qygridxdp=tunemap.outputs.Qygridxdp;
   
            Qxfmxdp(:,1:npd*npx) = [Qxgridxdp'-dqx; Qxgridxdp'-dqx; ...
                                    Qxgridxdp'+dqx; Qxgridxdp'+dqx];
            Qyfmxdp(:,1:npd*npx) = [Qygridxdp'-dqy; Qygridxdp'+dqy; ...
                                    Qygridxdp'+dqy; Qygridxdp'-dqy];
            figure;
            fill(Qxfmxdp,Qyfmxdp,Qdifxdp);hold on;
            xlabel('qx');ylabel('qy');
            plot_net(resorder,qxrange(1),qxrange(2),...
                          qyrange(1),qyrange(2));
            grid on;
            colormap('jet');
            shading flat;
            clim(caxrange);
            if (ratef)
                title(strcat(plottitle,{' Tune Diffusion Rate:(X,DP),Res order='},num2str(resorder)));
            else
                title(strcat(plottitle,{' Tune Diffusion:(X,DP), Res order='},num2str(resorder)));
            end
            colorbar;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qdifxdp not available for plottype %s \n', datetime, plottype);
        end

    case {'fmydp';'FMYDP'}
        if (ratef)
            Qdifydp   = tunemap.outputs.Qdifydpra;
        else
            Qdifydp   = tunemap.outputs.Qdifydp;
        end

        if(not(isempty(Qdifydp)))
            Qxgridydp=tunemap.outputs.Qxgridydp;
            Qygridydp=tunemap.outputs.Qygridydp;
   
            Qxfmydp(:,1:npd*npy) = [Qxgridydp'-dqx; Qxgridydp'-dqx; ...
                                    Qxgridydp'+dqx; Qxgridydp'+dqx];
            Qyfmydp(:,1:npd*npy) = [Qygridydp'-dqy; Qygridydp'+dqy; ...
                                    Qygridydp'+dqy; Qygridydp'-dqy];
            figure;
            fill(Qxfmydp,Qyfmydp,Qdifydp);hold on;
            xlabel('qx');ylabel('qy');
            plot_net(resorder,qxrange(1),qxrange(2),...
                          qyrange(1),qyrange(2));
            grid on;
            colormap('jet');
            shading flat;
            clim(caxrange);
            if (ratef)
                title(strcat(plottitle,{' Tune Diffusion Rate:(Y,DP),Res order='},num2str(resorder)));
            else
                title(strcat(plottitle,{' Tune Diffusion:(Y,DP),Res order='},num2str(resorder)));
            end
            colorbar;
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qdifydp not available for plottype %s \n', datetime, plottype);
        end

    case {'chro';'CHRO'}
         if not(isempty(dps)||isempty(Qxdpplot)||isempty(Qydpplot))
            figure; xlim([dpminplot dpmaxplot]*100);
            plot(dps*100,Qxdpplot,'-ok');xlabel('dp[%]');
            hold on;
            if (strcmpi(plotmode,'abs'))
                ylabel('Qx');
                yyaxis right; 
            else
                ylabel('dQ');
            end
            plot(dps*100,Qydpplot,'-or');
            if (strcmpi(plotmode,'abs'))
                ylabel('Qy');
                legend({'Qx';'Qy'});
            else
                legend({'dQx';'dQy'});
            end
            grid on;
            title(plottitle);
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qxdp or Qydp data not available in tunemap structure for plottype %s\n', datetime, plottype);
        end
        
    case {'chrotd';'CHROTD'}
        Qxdp=tunemap.outputs.Qxdp;
        Qydp=tunemap.outputs.Qydp;
        [~,dp0pos] = min(abs(dps));

        if not(isempty(Qxdp)||isempty(Qydp))
            figure;plot(Qxdp(1:dp0pos),Qydp(1:dp0pos),'ob');hold on;
            plot(Qxdp(dp0pos:length(dps)),Qydp(dp0pos:length(dps)),'or');
            plot_net(resorder,qxrange(1),qxrange(2),...
                          qyrange(1),qyrange(2));
                          
            xlabel('qx');ylabel('qy');
            title(strcat(plottitle,{' Res order='},num2str(resorder)));
            nhandles=nhandles+1;
            phandles{nhandles}=gcf;
        else
            fprintf('%s Error in plotTuneMap: Qxx, Aqxy, Qyx or Qyx not available for plottype %s \n', datetime, plottype);
        end
    otherwise
        fprintf('%s Error in plotTuneMap: unknown plottype: %s \n',datetime,plottype);       
        return
end
%% Auxiliary functions
function plotslice(ax,ay,Q,plane,z0)
%
% plots vertical or horizontal slices of a two-dimensional tune 
%
% inputs
% ax : (1Xn) matrix of horizontal coordinates
% ay : (1Xm) matrix of vertical coordinates
% Q  : (n*mX1) matrix of values to be plotted, grouped by rows,

n=numel(ax);
m=numel(ay);
[amplx_m, amply_m]=meshgrid(ax,ay);
Q_m = reshape(Q,m,n);
legstr=cell(1,numel(z0));

% handles the case of z0='all'
if (ischar(z0))
    switch plane
        case {'x';'X'}
            z0=ay;
        case {'y';'N'}
            z0=ax;
    end
end

for i=1:numel(z0)
    switch plane
        case {'x';'X'}  
            [~,idx] = min(abs(ay-z0(i)));
            legstr{i}=num2str(ay(idx));
            if (i==1)
                figure;plot(amplx_m(idx,:),Q_m(idx,:),'o-');
                hold on;
            else
                plot(amplx_m(idx,:),Q_m(idx,:),'o-');
            end

        case {'y','Y'}
            [~,idx] = min(abs(ax-z0(i)));
            legstr{i}=num2str(ax(idx));
            if (i==1)
                figure;plot(amply_m(:,idx),Q_m(:,idx),'o-');hold on;
            else
                plot(amply_m(:,idx),Q_m(:,idx),'o-');
            end
    end
end
legend(legstr);
hold off;



function plotDA(varargin)
% plots the Dynamic aperture
%   
%% Inputs
% Mandatory argument
% DAS : structure generated by calcDA
% 
% Optional arguments
% xminplot : minimum horizontal limit for DA plot [m]
% xmaxplot : maximum horizontal limit for DA plot [m]
% ymaxplot : maximum vertical limits for DA plots are (0,ymaxplot) [m]
%            If not given, the two values above are taken from
%            DAS.outputs.DAoptions
%
% dpminplot : minimum energy deviation for xdp and ydp plots
% dpmaxplot : maximum energy deviation for xdp and ydp plots
%
% verbose : defines level of verbose output, default=0, i.e. no output
%
%% Usage examples
% plotDA(DAS);
% plotDA(DAS,'xminplot', -0.015, 'xmaxplot', 0.012, 'ymaxplot', 0.004);

%% History
% PFT 2024/05/01
% PFT 2024/05/05 bug fix for grid mode plotting
% PFT 2024/05/14 added xp and yp plots
% PFT 2024/05/25 updated handling of plot limit parameters
% PFT 2024/05/26 updated handling of verbose output

%% Input argument parsing
%
DAS       = getargs(varargin,[]);
verbosef  = getoption(varargin,'verbose',0);
xminplot  = getoption(varargin,'xminplot',-DAS.outputs.DAoptions.XmaxDA);
xmaxplot  = getoption(varargin,'xmaxplot', DAS.outputs.DAoptions.XmaxDA);
ymaxplot  = getoption(varargin,'ymaxplot',DAS.outputs.DAoptions.YmaxDA);
dpminplot = getoption(varargin,'dpminplot',DAS.outputs.DAoptions.dpmin);
dpmaxplot = getoption(varargin,'dpminplot',DAS.outputs.DAoptions.dpmax);

DAmode = DAS.outputs.DAoptions.DAmode;
dp     = DAS.outputs.DAoptions.dp;
DAV    = DAS.outputs.DAV;
npdax  = DAS.outputs.DAoptions.npdax;
npday  = DAS.outputs.DAoptions.npday;

if (not(isfield(DAS.inputs,'mode')))
    DAS.inputs.mode='xy';
end

switch DAS.inputs.mode
    case {'xy';'XY'}    
        switch DAmode
            case 'border'
                figure;plot(DAV(:,1)*1000,DAV(:,2)*1000,'-ob');
                xlabel('X [mm]'); ylabel('Y [mm]');grid;
                xlim([xminplot xmaxplot]*1000);ylim([0 ymaxplot]*1000);
                title(sprintf('dp = %3.1f %%', dp*100));
                grid on;

            case 'grid'
                DAM = zeros(npday+1,2*npdax+1);
                k= 1;
                for i=0:npday
                    for j= 1:2*npdax+1
                        DAM(npday+1-i,j)=DAV(k);
                        k=k+1;
                    end
                end
                DAM=DAM*255;
                map=[0 0.75 0; 1 1 1];
                figure;image([xminplot xmaxplot]*1000,[Ymaxplot*1000,0],DAM);
                ax=gca;
                ax.YDir='normal';
                colormap(map);
                xlabel('X[mm]');
                ylabel('Y[mm]');    
                xlim([xminplot xmaxplot]*1000);ylim([0 ymaxplot]*1000);grid;  
                title(sprintf('dp = %3.1f %%', dp*100));
                grid on;
        end

    case {'xydp';'XYDP'}
        dps=DAS.outputs.dps;
        DAXp=DAS.outputs.DAXp;

        DAYp=DAS.outputs.DAYp;
        DAXm=DAS.outputs.DAXm;

        figure;plot(dps*100,DAXp*1000,'-ob');hold;plot(dps*100,DAXm*1000,'-or');
        xlabel('dp[%]');ylabel('X[mm]');
        xlim([dpminplot dpmaxplot]*100);ylim([xminplot xmaxplot]*1000);
        grid on;

        figure;plot(dps*100,DAYp*1000,'-ob');
        xlabel('dp[%]');ylabel('Y[mm]');
        xlim([dpminplot dpmaxplot]*100);ylim([0 ymaxplot]*1000);
        grid on;
end



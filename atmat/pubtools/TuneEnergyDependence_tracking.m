function [nux, nuy, MDnux, MDnuy] = TuneEnergyDependence_tracking(RING,nturn,ndp,dppmin,dppMAX)

% --------------------------------------------------------
% calculate tune shift with energy using tracking and NAFF
%
% INPUT:
% RING        = an AT2 compliant lattice ring (def: RING)
% nturn       = the number of turns used to compute tunes and diffusion
% ndp         = number of energy bins
% dppmin,MAX  = min/MAX value of energy (should be symmetric)
% 
% OUTPUT:
% nux,y       = the tunes as a function of energy shift
% MDnux,y     = maximal displacement of the tunes from the nominal (dpp=0)
%               value
% --------------------------------------------------------

c = parcluster;
% next line is optional at MAX IV
c.AdditionalProperties.AccountName = 'any-virtual-account-name';
% 6 hour walltime
c.AdditionalProperties.WallTime = '06:00:00';
% hyperthreading enabled
c.NumThreads = 6;
c.saveProfile;
%%%%parpool('aurora R2022a',56) %%% IMPORTANT REFERENCE ---
parpool('local',12);
pp = gcp; 

if nargin<2
nturn = 2048; 
ndp   = 17;
dppmin = -0.04;
dppMAX = +0.04;
end

RINGt = atdisable_6d(RING); % remove cavity/radiation effects and control energy "manually"
dpp = linspace(dppmin,dppMAX,ndp);
dx = 1e-4; dy = 1e-4; % small dx,y perturbation to allow the NAFF calculation

parfor i = 1:numel(dpp)
    Ri    = [dx 0 dy 0 dpp(i) 0]';
    Ro{i} = ringpass(RINGt, Ri, nturn);
end

nfreq=1; % only one freq to process
for i = 1:numel(dpp)
    vv(1:6,1,1:nturn) = Ro{i}(1:6,1:nturn);
    v  = vv;
    xf = squeeze(v(:,1,:));
    mxhalf1 = mean(xf(:,1:nturn)');
    xhalf   = xf(:,1:nturn)-mxhalf1'*ones(1,nturn); % subtracting mean value
    [freq amplitude phase] = calcnaff(xhalf(1,:), xhalf(2,:));
    [a b] = max(abs(amplitude));
    %nf = length(freq); if (nf < nfreq), freq((nf+1):nfreq) = 0.0; end;
    nux(i) = abs(freq(b)/2/pi);

    [freq amplitude phase] = calcnaff(xhalf(3,:), xhalf(4,:));
    [a b] = max(abs(amplitude));
    %nf = length(freq); if (nf < nfreq), freq((nf+1):nfreq) = 0.0; end;
    nuy(i) = abs(freq(b)/2/pi);
end
pdm   = numel(nux)/2+mod(numel(nux),2)/2;
nux0 = nux(pdm);
nuy0 = nuy(pdm);   

if mod(numel(nux),2)==0
   nux0 = mean([nux(pdm) nux(pdm+1)]);
   nuy0 = mean([nuy(pdm) nuy(pdm+1)]);   
end

MDnux = max(abs(nux-nux0)); % maximal tune excursion from 
MDnuy = max(abs(nuy-nuy0)); % the nominal point at dpp=0
                                % can be used for ADTS optimsations 

delete(pp)

end
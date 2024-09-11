function smm = atsummary_fast(ring,isdipole)


%    e_mass=510999.06450473;	% eV
%    e_radius= 2.81794091754E-15;                  % m
%    hbar=197.327054064; % MeV*fm
%    cspeed = 2.99792458e8;                   % m/s
%    Cq=55/32/sqrt(3)*hbar/e_mass*1E-9; % m
     Cq = 3.831939211875470e-13;
%    Cgamma=4.0E27*pi*e_radius/3/e_mass^3;   % [m/GeV^3]
%    e0=3; % energy eV
     gamma=5870.85223514;

try
    [ringdata,TD] = atlinopt4(ring,1:length(ring)+1,'get_chrom','coupled',false);

    smm.circ = findspos(ring,length(ring)+1);
    smm.tunes = [TD(length(TD)).mu(1)/(2*pi), TD(length(TD)).mu(2)/(2*pi)];
    smm.Qx_ring = smm.tunes(1)*20;
    smm.Qy_ring = smm.tunes(2)*20;
    smm.chromaticity = ringdata.chromaticity;

    % For calculating the synchrotron integrals
    [~,I2d,I3d,I4d,I5d,~,~] = DipoleRadiation_fast(ring,TD,isdipole);
    %smm.integrals=[I1d,I2d,I3d,I4d,I5d,I6];

    alphac=mcf(ring,0.0);
    %eloss=1.0e9*Cgamma/2/pi*e0.^4*smm.integrals(2); % eV    
    %smm.compactionFactor = alphac;

    % Damping numbers
    % Use Robinson's Theorem
    %smm.damping(1) = 1 - smm.integrals(4)/smm.integrals(2);
    smm.damping(1) = 1 - I4d/I2d;
    %smm.damping(2) = 1;
    smm.damping(3) = 2 + I4d/I2d;

    %smm.radiation = 1.0e-9*eloss;   % GeV

    smm.naturalEnergySpread = gamma*sqrt(Cq*I3d/(2*I2d + I4d));
    smm.naturalEmittance    = Cq*gamma.^2*I5d/(I2d-I4d);
    smm.alphac = alphac;
    smm.etax = TD(1).Dispersion(1);
    smm.beta0 = TD(1).beta;

%[ring,newargs]=getargs(varargin,THERING,'check',@iscell);
%[varargout{1:nargout}]=wrapper6d(ring,@doit,newargs{:});
catch 
    smm.circ = nan;
    smm.tunes(1)= nan;
    smm.tunes(2)= nan;
    smm.Qx_ring = nan;
    smm.Qy_ring = nan;
    smm.chromaticity(1) = nan;
    smm.chromaticity(2) = nan;
    smm.damping(1) = nan;
    smm.damping(2) = nan;
    smm.naturalEnergySpread = nan;
    smm.naturalEmittance    = nan;
    smm.alphac = nan;
    smm.etax = nan;
    smm.beta0(1) = nan;
    smm.beta0(2) = nan;
end


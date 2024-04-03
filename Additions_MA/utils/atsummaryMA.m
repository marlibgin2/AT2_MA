function Summary = atsummaryMA
%ATSUMMARY - Prints out the paramters of the current AT lattice
%  The parameters that come after the Synchrotron Integrals are
%  parameters that depend on the Integrals themselves. The equations to
%  calculate them were taken from [1].
%
%  [1] Alexander Wu Chao and Maury Tigner, Handbook of Accelerator Physics
%  and Engineering (World Scientific, Singapore, 1998), pp. 183-187. (or
%  187-190 in ed. 2)
%
%  See also ringpara

%  Written by Eugene Tan
%  Revised by Laurent S. Nadolski


global THERING

% Structure to store info
Summary.e0 = THERING{1}.Energy/1e9; %getenergy('Model');
Summary.circumference = findspos(THERING, length(THERING)+1);
Summary.revTime = Summary.circumference / 2.99792458e8;
Summary.revFreq = 2.99792458e8 / Summary.circumference;
Summary.gamma = Summary.e0 / 0.51099906e-3;
Summary.beta = sqrt(1 - 1/Summary.gamma);
[TD, Summary.tunes, Summary.chromaticity] = twissring(THERING, 0, 1:length(THERING)+1, 'chrom', 1e-8);
Summary.compactionFactor = mcf(THERING);

% For calculating the synchrotron integrals
temp  = cat(2,TD.Dispersion);
D_x   = temp(1,:)';
D_x_  = temp(2,:)';
beta  = cat(1, TD.beta);
alpha = cat(1, TD.alpha);
gamma = (1 + alpha.^2) ./ beta;
circ  = TD(length(THERING)+1).SPos;

% Synchrotron integral calculation
Summary.integrals = [0.0 0.0 0.0 0.0 0.0 0.0];

for i = 1:length(THERING)
    if isfield(THERING{i}, 'BendingAngle') && isfield(THERING{i}, 'EntranceAngle')
        rho = THERING{i}.Length/THERING{i}.BendingAngle;
        dispersion = 0.5*(D_x(i)+D_x(i+1));
        Summary.integrals(1) = Summary.integrals(1) + dispersion*THERING{i}.Length/rho;
        Summary.integrals(2) = Summary.integrals(2) + THERING{i}.Length/(rho^2);
        Summary.integrals(3) = Summary.integrals(3) + THERING{i}.Length/(rho^3);
        % For general wedge magnets
        Summary.integrals(4) = Summary.integrals(4) + ...
            D_x(i)*tan(THERING{i}.EntranceAngle)/rho^2 + ...
            (1 + 2*rho^2*THERING{i}.PolynomB(2))*(D_x(i)+D_x(i+1))*THERING{i}.Length/(2*rho^3) + ...
            D_x(i+1)*tan(THERING{i}.ExitAngle)/rho^2;
        %         Summary.integrals(4) = Summary.integrals(4) + 2*0.5*(D_x(i)+D_x(i+1))*THERING{i}.Length/rho^3;
        H1 = beta(i,1)*D_x_(i)*D_x_(i)+2*alpha(i)*D_x(i)*D_x_(i)+gamma(i)*D_x(i)*D_x(i);
        H0 = beta(i+1,1)*D_x_(i+1)*D_x_(i+1)+2*alpha(i+1)*D_x(i+1)*D_x_(i+1)+gamma(i+1)*D_x(i+1)*D_x(i+1);
        Summary.integrals(5) = Summary.integrals(5) + THERING{i}.Length*(H1+H0)*0.5/(rho^3);
        %         if H1+H0 < 0
        %             fprintf('%f %i %s\n', H1+H0, i, THERING{i}.FamName)
        %         end
        Summary.integrals(6) = Summary.integrals(6) + THERING{i}.PolynomB(2)^2*dispersion^2*THERING{i}.Length;
    end
end

% Damping numbers
% Use Robinson's Theorem
Summary.damping(1) = 1 - Summary.integrals(4)/Summary.integrals(2);
Summary.damping(2) = 1;
Summary.damping(3) = 2 + Summary.integrals(4)/Summary.integrals(2);

Summary.radiation = 8.846e-5*Summary.e0.^4*Summary.integrals(2)/(2*pi);
Summary.naturalEnergySpread = sqrt(3.8319e-13*Summary.gamma.^2*Summary.integrals(3)/(2*Summary.integrals(2) + Summary.integrals(4)));
Summary.naturalEmittance = 3.8319e-13*(Summary.e0*1e3/0.510999).^2*Summary.integrals(5)/(Summary.damping(1)*Summary.integrals(2));

% Damping times
Summary.radiationDamping(1) = 1/(2113.1*Summary.e0.^3*Summary.integrals(2)*Summary.damping(1)/circ);
Summary.radiationDamping(2) = 1/(2113.1*Summary.e0.^3*Summary.integrals(2)*Summary.damping(2)/circ);
Summary.radiationDamping(3) = 1/(2113.1*Summary.e0.^3*Summary.integrals(2)*Summary.damping(3)/circ);

% Slip factor
Summary.etac = Summary.gamma^(-2) - Summary.compactionFactor;

cavind = findcells(THERING,'HarmNumber');
if ~isempty(cavind)
    freq    = THERING{cavind(1)}.Frequency;
    v_cav   = getcellstruct(THERING,'Voltage',cavind);
    I       = findcells(THERING,'TimeLag');             % TimeLag is used instead of PhaseLag in the CavityPass passmethod, but might not be defined
    if ~isempty(I)
        cT  = zeros(size(v_cav));
        [~, I, ~] = union(I,cavind);
        cT(I)  = getcellstruct(THERING,'TimeLag',cavind(I));
    else
        cT  = zeros(size(v_cav));
    end
%     phi = getcellstruct(THERING,'PhaseLag',cavind);
    phi     = 2*pi * cT ./ (299792458 ./ THERING{cavind(1)}.Frequency);
    v_cav   = abs(sum(v_cav.*sin(phi) + sqrt(-1)*v_cav.*cos(phi)));
else
    % Default
    freq = 352.202e6;
    v_cav = 3e6;
end
Summary.harmon = Summary.circumference/(2.99792458e8/freq); % AsSuming 499.654MHz RF
Summary.overvoltage = v_cav/(Summary.radiation*1e9); % AsSuming 3e6 volt cavities.
% Assuming the harmon and overvoltage above.
% references:  H. Winick, "Synchrotron Radiation Sources: A Primer",
% World Scientific Publishing, Singapore, pp92-95. (1995)
% Wiedemann, pp290,350. Chao, pp189.
Summary.syncphase = pi - asin(1/Summary.overvoltage);
Summary.energyacceptance = sqrt(v_cav*sin(Summary.syncphase)*2*(sqrt(Summary.overvoltage^2-1) - acos(1/Summary.overvoltage))/(pi*Summary.harmon*abs(Summary.etac)*Summary.e0*1e9));
Summary.synctune = sqrt((Summary.etac*Summary.harmon*v_cav*cos(Summary.syncphase))/(2*pi*Summary.e0*1e9));
Summary.bunchlength = Summary.beta*299792458*abs(Summary.etac)*Summary.naturalEnergySpread/(Summary.synctune*Summary.revFreq*2*pi);

if nargout == 0
    fprintf('\n');
    %fprintf('   ******** Summary for ''%s'' ********\n', GLOBVAL.LatticeFile);
    fprintf('   ******** AT Lattice Summary ********\n');
    fprintf('   Energy: \t\t\t%4.5f [GeV]\n', Summary.e0);
    fprintf('   Gamma: \t\t\t%4.5f \n', Summary.gamma);
    fprintf('   Circumference: \t\t%4.5f [m]\n', Summary.circumference);
    fprintf('   Revolution time: \t\t%4.5f [ns] (%4.5f [MHz]) \n', Summary.revTime*1e9,Summary.revFreq*1e-6);
    fprintf('   Betatron tune H: \t\t%4.5f (%4.5f [kHz])\n', Summary.tunes(1),Summary.tunes(1)/Summary.revTime*1e-3);
    fprintf('                 V: \t\t%4.5f (%4.5f [kHz])\n', Summary.tunes(2),Summary.tunes(2)/Summary.revTime*1e-3);
    fprintf('   Momentum Compaction Factor: \t%4.5f\n', Summary.compactionFactor);
    fprintf('   Chromaticity H: \t\t%+4.5f\n', Summary.chromaticity(1));
    fprintf('                V: \t\t%+4.5f\n', Summary.chromaticity(2));
    fprintf('   Synchrotron Integral 1: \t%4.5f [m]\n', Summary.integrals(1));
    fprintf('                        2: \t%4.5f [m^-1]\n', Summary.integrals(2));
    fprintf('                        3: \t%4.5f [m^-2]\n', Summary.integrals(3));
    fprintf('                        4: \t%4.5f [m^-1]\n', Summary.integrals(4));
    fprintf('                        5: \t%4.5f [m^-1]\n', Summary.integrals(5));
    fprintf('                        6: \t%4.5f [m^-1]\n', Summary.integrals(6));
    fprintf('   Damping Partition H: \t%4.5f\n', Summary.damping(1));
    fprintf('                     V: \t%4.5f\n', Summary.damping(2));
    fprintf('                     E: \t%4.5f\n', Summary.damping(3));
    fprintf('   Radiation Loss: \t\t%4.5f [keV]\n', Summary.radiation*1e6);
    fprintf('   Natural Energy Spread: \t%4.5e\n', Summary.naturalEnergySpread);
    fprintf('   Natural Emittance: \t\t%4.5e [mrad]\n', Summary.naturalEmittance);
    fprintf('   Radiation Damping H: \t%4.5f [ms]\n', Summary.radiationDamping(1)*1e3);
    fprintf('                     V: \t%4.5f [ms]\n', Summary.radiationDamping(2)*1e3);
    fprintf('                     E: \t%4.5f [ms]\n', Summary.radiationDamping(3)*1e3);
    fprintf('   Slip factor : \t%4.5f\n', Summary.etac);
    fprintf('\n');
    fprintf('   Assuming cavities Voltage: %4.5f [kV]\n', v_cav/1e3);
    fprintf('                   Frequency: %4.5f [MHz]\n', freq/1e6);
    fprintf('             Harmonic Number: %4.5f\n', Summary.harmon);
    fprintf('   Overvoltage factor: %4.5f\n', Summary.overvoltage);
    fprintf('   Synchronous Phase:  %4.5f [rad] (%4.5f [deg])\n', Summary.syncphase, Summary.syncphase*180/pi);
    fprintf('   Linear Energy Acceptance:  %4.5f %%\n', Summary.energyacceptance*100);
    fprintf('   Synchrotron Tune:   %4.5f (%4.5f kHz or %4.2f turns) \n', Summary.synctune, Summary.synctune/Summary.revTime*1e-3, 1/Summary.synctune);
    fprintf('   Bunch Length:       %4.5f [mm]\n', Summary.bunchlength*1e3);
end
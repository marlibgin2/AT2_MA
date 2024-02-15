function rbpars = RBParam(dTheta, varargin)
%RBParam 
%   Returns bem parameters for ring composed only of unit cells
%   with two reverse bends.
%   Input: dTheta 
%          dTheta = Reverse bend angle  [rad]
%          varargin{1} reverse bend quadrupolar strength - default is same
%          as in lattice without reverse bend [m**-2]
%
%   RUCRB is a global cell array containing the lattice
%
%   Ouput: [dTheta*1000 Krb Circ dCirc dfRF 0.0 yOffset unpkick RbShift ...
%           emit Remit Jx alphac sigma_E Qx Qy Ksix Ksiy betax0 betay0 disp0];
%   dTheta = Reverse bend angle [mrad]
%   Krb    = Reverse Bend Quadrupole strength [m**-2]
%   Circ   = Circumference[m]
%   dCirC  = change in circumference compare to zero RB kick case [mm]
%   dfRF   = change in RF frequency [Hz]
%   TotAngle = total bend angle [deg]
%   yOffset  = vertical orbvoit offset at the position of the Reverse Bends [mm]
%   unpkick  = kick from Reverse bend due to the vertical offset above [mrad]
%   RbShift  = additional horiuzonatl shift of Reverse bend in order to
%   achieve requierd total RB kick [mm]
%   emit     = emittance [pm rad]
%   Remit    = Ration of emittance to unperturbed case [%]
%   Jx       = Horizontal dampinmg partition
%   alphac   = momentum compaction [0.001]
%   sigma_E  = energy spread [0.001]
%   Qx, Qy   = Betatron Tunes
%   Ksix Ksiy = Chormaticities
%   betax0 betay0 disp0 = betatron and dispersion fucntions at center of dipole

global RUCRB UCRB


frf=99.931E6;
RbLength = 0.15;

if (nargin>1)
    Krb = varargin{1};
    max4_20121107_430e_studies('Yes',dTheta,Krb);
else
    max4_20121107_430e_studies('Yes',dTheta);
    Krb = 4.030076;
end

[x2d y2d a2d] = Survey2D(UCRB,1.5*pi/180);
yOffset = 6.509501134666746-y2d(round(length(y2d)/2))*1000;
pars=ringpara(RUCRB);
Circ = pars.R*2*pi;
dCirc = (Circ-336.000)*1000;
dfRF  = dCirc/336*frf/1000;
unpkick = RbLength*Krb*yOffset;
RbShift = (dTheta*1000 - unpkick)/RbLength/Krb;
emit = pars.emittx*1e12;
Remit = emit/336.9050205153787*100;
Qx = pars.tunes(1);
Qy = pars.tunes(2);
Jx  = pars.dampingJ(1);
Ksix = pars.chroms(1);
Ksiy = pars.chroms(2);
alphac = pars.alphac*1000;
sigma_E = pars.sigma_E*1000;

tw=twissring(UCRB,0,[1],'chroma');
betax0 = tw.beta(1);
betay0 = tw.beta(2);
disp0  = tw.Dispersion(1);

rbpars = [dTheta*1000 Krb Circ dCirc dfRF 0.0 yOffset unpkick RbShift ...
           emit Remit Jx alphac sigma_E Qx Qy Ksix Ksiy betax0 betay0 disp0];

end


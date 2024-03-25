function [c,ceq] = stabcon(x,LatticeOptData)
%
% inputs: x - lattice devision variables
%         Lattice OptData: struture with Lattice data - generated by
%                          script (max4_UpgradeStudies).
% outputs: c(1)=|Trace(Mx)|-2. 
%          c(2)=|Trace(Mx)|-2
%          c(3)= Maxbetax - BetaXMAX;
%          c(4)= Maxbetay - BetaYMAX;
%          c(5)= abs(FractuneGoal_X-FracTuneX)-DFracTuneX
%          c(6)= abs(FractuneGoal_X-FracTuneX)-DFracTuneX
%          c(7)= abs(etax0) - EtaX0MAX
%          c(8)= BetaY0 - BetaX0 + DBetaXY
%          c(9)= BetaX0 - BetaX0MAX
%          c(10)=BetaY0 - BetaY0MAX
%          c(11)= Jx - 3.0
%
%   constraint outpts must be < 0 for a lattice to be accepted (feasible)
%   Stable latices require c(1),c(2)<0 
%
%   if BetaXMax, BetaYMAX, EtaX0, FracTuneGoal_X, FracTuneGoal_Y or DBetaXY 
%       are NaN, the corresponding conditions are not required, i.e. the output values
%       are set to -1.
%
%   Note: for conditions 5 to 8 to be properly evaluated, BetXMAX and
%   BetaYMax must not be nan, i.e. they must be evaluated also !
% 
%

HACHRO      = LatticeOptData.HACHRO;
ACHRO       = LatticeOptData.ACHRO;
BetaXMAX    = LatticeOptData.BetaXMAX;
BetaYMAX    = LatticeOptData.BetaYMAX;
EtaX0MAX    = LatticeOptData.EtaX0MAX;
GoalFracTuneX = LatticeOptData.GoalFracTuneX;
GoalFracTuneY = LatticeOptData.GoalFracTuneY;
DeltaFracTuneX = LatticeOptData.DeltaFracTuneX;
DeltaFracTuneY = LatticeOptData.DeltaFracTuneY;
isdipole = LatticeOptData.isdipole;

if (isfield(LatticeOptData,'DBetaXY'))
    DBetaXY        = LatticeOptData.DBetaXY;
else
    DBetaXY = nan;
end

if (isfield(LatticeOptData,'BetaX0Max'))
    BetaX0Max        = LatticeOptData.BetaX0Max;
else
    BetaX0Max = nan;
end

if (isfield(LatticeOptData,'BetaY0Max'))
    BetaY0Max        = LatticeOptData.BetaY0Max;
else
    BetaY0Max = nan;
end

HACHRO = setDVs(1,HACHRO,LatticeOptData, x);
ACHRO  = setDVs(2, ACHRO,LatticeOptData, x);
%
% Evaluates constraints
%
props=LatticeOptData.props;
ceq=[];
RIN=LatticeOptData.RIN;
NE=LatticeOptData.NE;
scaling=LatticeOptData.scaling;

ROUT=atpass(HACHRO,RIN,1,1,NE+1,cell(0),cell(0),1,1,props);

TMAT3 = reshape(ROUT(1:4,:),4,9,[]);
M44 = (TMAT3(:,1:4,end)-TMAT3(:,5:8,end))./scaling;
TRx = 2*(M44(1,1)*M44(2,2)+M44(1,2)*M44(2,1));
TRy = 2*(M44(3,3)*M44(4,4)+M44(3,4)*M44(4,3));
c(1) = abs(TRx)-2.0;
c(2) = abs(TRy)-2.0;

if ( (c(1)<0)&&(c(2)<0) )
    if (isnan(BetaXMAX)&&isnan(BetaYMAX))
        c(3)=-1;
        c(4)=-1;
    else
        [RD,TD] = atlinopt4(ACHRO,1:length(ACHRO)+1,'coupled',false);
        [~,I2d,~,I4d,~,~,~] = DipoleRadiation_fast(ACHRO,TD,isdipole);
        Jx = 1 - I4d/I2d;
        betas=cat(1,TD.beta);
        Maxbetax = max(betas(:,1));
        Maxbetay = max(betas(:,2));
        BetaX0 = betas(1,1);
        BetaY0 = betas(1,2);

        c(3)= Maxbetax - BetaXMAX;
        c(4)= Maxbetay - BetaYMAX;
        c(11) = Jx-3.0;

        fptx=(RD.tune(1)*20-fix(RD.tune(1)*20));
        fpty=(RD.tune(2)*20-fix(RD.tune(2)*20));
        if (isnan(GoalFracTuneX)&&isnan(GoalFracTuneY))
            c(5)=-1;
            c(6)=-1;
        else
            c(5)= abs(GoalFracTuneX-fptx)-DeltaFracTuneX;
            c(6)= abs(GoalFracTuneY-fpty)-DeltaFracTuneY;
        end
    end

    if (isnan(EtaX0MAX))
        c(7) = -1;
    else
        Etax0 = abs(TD(1).Dispersion(1));
        c(7) = Etax0 - EtaX0MAX;
    end

    if (isnan(DBetaXY))
        c(8) = -1;
    else
        c(8) = BetaY0 - BetaX0 + DBetaXY;
    end

    if (isnan(BetaX0Max))
        c(9)=-1;
    else
        c(9) = BetaX0 - BetaX0Max;
    end

    if (isnan(BetaY0Max))
        c(10)=-1;
    else
        c(10) = BetaY0 - BetaY0Max;
    end

else
    c(3)= 1;
    c(4)= 1;
    c(5)= 1;
    c(6)= 1;
    c(7)= 1;
    c(8)= 1;
    c(9)= 1;
    c(10)=1;
    c(11)=1;
end
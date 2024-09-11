function [NewLAT, its, penalty, ftunes] = fittuneRS(LAT,tunes, fam1, fam2, varargin)
% 
%   This is a wrapper function that iterates calls to fittune until 
%   a given tolerance or a maximum numer of iterations is achieved
%
%% Inputs
% Mandatory arguments
%
% LAT:   AT2.0 lattice structure
% tunes: (1x2) array of target tunes
% fam1, fam2: quadrupole families to be used in the fit
%
% Optional arguments
% maxits : maximum number of iterations, default = 10
% Tol    : tolerance for tune deviation, defaule = 1.0E-4
%          penalty function is defined as
%          (sqrt(targettunex-tunex)**2 + (targettuney-tuney)**2)
% frac   : fraction of tune correction to implement at each iteration,
%          default=1.0
% UseIntegarPart: fits the tune using the integer part, default=true;
% verbose: defines level of verbose output (default = 0)
%
%
% %% Outputs
% NewLat : new lattice with corrected tunes
% its : number of actual iterations
% penalty : achieved penalty function 
% ftunes : final tunes
%
%% Usage examples
% [RINGc,its,penalty,~]=fittuneRS(RINGe,[55.398 15.199],'Q1','Q2','maxits',10,'Tol',1E-5,'frac',0.9,'UseIntegerPart')
%

%% History
% PFT: sometime in 2023, first version
% PFT 2024/05/24: general review,improved documentation, added frac as
%                 optional argument
% PFT 2024/05/25: changed the way to handle verbose output level
%
%% Input argument parsing
maxits           = getoption(varargin,'maxits',10);
Tol              = getoption(varargin,'Tol',1.0E-4);
frac             = getoption(varargin,'frac',1.0);
verboselevel     = getoption(varargin,'verbose',0);
useintegerpartf  = getoption(varargin,'UseIntegerPart',true);


%% Preamble
Ifam1=find(atgetcells(LAT, 'FamName', fam1),1);
Ifam2=find(atgetcells(LAT, 'FamName', fam2),1);
if (isempty(Ifam1)||isempty(Ifam2))
    if (isempty(Ifam1))
         fprintf('%s Error in fittuneRS qpole family %s not in lattice \n', datetime, fam1);
    end
    if (isempty(Ifam2))
        fprintf('%s Error in fittuneRS qpole family %s not in lattice \n', datetime, fam2);
    end
    NewLAT=LAT;
    its=nan;
    penalty=nan;
    ftunes=[nan nan];
    return
end

if useintegerpartf
   gettune = @getinttune;
else
   if any(tunes>=1)
       if (verboselevel>0)
            fprintf('%s FittruneRS - The integer part of the tunes is ignored unless you use the ''UseIntegerPart'' flag \n', datetime);
            tunes=tunes-floor(tunes);
       end
       gettune = @getfractune;
   end
end
tunes_now=gettune(LAT);

penalty=sqrt(sum((tunes-tunes_now).^2));

%% Fits the tunes
NewLAT=LAT;
for i=1:maxits
    if (penalty<Tol); break
    else
       NewLAT=fittuneR(NewLAT, tunes, fam1, fam2, 'UseIntegerPart', useintegerpartf, 'frac',frac); 
       tunes_now=gettune(NewLAT);
       %switch useintegerpartf
       %    case 'Y'
       %         NewLAT=fittuneR(NewLAT, tunes, fam1, fam2, 'UseIntegerPart', 'frac',frac);
       %         tunes_now=getinttune(NewLAT);
       %    otherwise
       %         NewLAT=fittuneR(NewLAT, tunes, fam1, fam2,'frac',frac);
       %         tunes_now=getfractune(NewLAT);
       %end
       penalty=sqrt(sum((tunes-tunes_now).^2));
    end
end
its=i;
ftunes = tunes_now;
end
%% Auxiliary functions
function tun=getfractune(ring,varargin)
            tun3=tunechrom(ring,varargin{:});
            tun=tun3(1:2);
end
        
function tun=getinttune(ring)
            [~,TD] = atlinopt4(ring,1:length(ring)+1,'coupled',false);
            tun= TD(end).mu(1:2)/2/pi;
end


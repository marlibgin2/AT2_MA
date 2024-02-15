function [NewLAT, its, penalty, ftunes] = fittuneRS(LAT,tunes, fam1, fam2, maxits, Tol,useintegerf)
% Iterates calls to fittuneR
%   Detailed explanation goes here

switch useintegerf
    case 'Y'
    tunes_now=getinttune(LAT);
    otherwise
    tunes_now=getfractune(LAT);
end
penalty=sum((tunes-tunes_now).^2);

NewLAT=LAT;
for i=1:maxits
    if (penalty<Tol); break
    else
       switch useintegerf
           case 'Y'
                NewLAT=fittuneR(NewLAT, tunes, fam1, fam2, 'UseIntegerPart');
                tunes_now=getinttune(NewLAT);
           otherwise
                NewLAT=fittuneR(NewLAT, tunes, fam1, fam2);
                tunes_now=getfractune(NewLAT);
       end
       penalty=sum((tunes-tunes_now).^2);
    end
end
its=i;
ftunes = tunes_now;
end

function tun=getfractune(ring,varargin)
            tun3=tunechrom(ring,varargin{:});
            tun=tun3(1:2);
end
        
  function tun=getinttune(ring)
            [~,TD] = atlinopt4(ring,1:length(ring)+1,'coupled',false);
            tun= TD(end).mu(1:2)/2/pi;
        end
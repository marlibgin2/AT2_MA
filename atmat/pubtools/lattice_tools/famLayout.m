function famLayout=famLayout(ACHRO)
%
% Returns a list of magnet family names in the sequence they appear in the
% input lattice cell array
%% Inputs
% ACHRO : AT2 lattice cell array
%
%% Outputs
% famLayout (NX3): cell array with two coluns. The first 
%              contains strings with names of magnet families in
%              the order they appear in the lattice. The 
%              second contains the corresponding multipolar strengths 
%              of those element - the PolynomB array field of the element.
%              The third contains the bend angle of the element in degrees
% 
%% Usage examples
% fLO=famLayout(ACHRO);
%

%% History
% PFT 2024/07/17 : firs version
% PFT 2024/07/24 : added girder markers, correctors and monitors 
%                  if available
% PFT 2024/08/03 : fixed bug update counter for correctors
%
famLayout={};
j=0;
for i=1:numel(ACHRO)
    element=ACHRO{i};
    if (isfield(element,'PassMethod'))
        PM=element.PassMethod;
        if (strcmp(PM,'BndMPoleSymplectic4Pass')||strcmp(PM,'StrMPoleSymplectic4Pass'))
            if (isfield(element,'BendingAngle'))
                bendangle=element.BendingAngle;
            else
                bendangle=0.0;
            end
            if (j>0)
                if(not(strcmp(element.FamName,famLayout{j,1})))
                    famLayout=[famLayout;{element.FamName},element.PolynomB,bendangle*180/pi];
                    j=j+1;
                end
            else
                famLayout=[famLayout;{element.FamName},element.PolynomB,bendangle*180/pi];
                j=j+1;
            end
        end
        if (strcmp(PM,'CorrectorPass'))
            famLayout=[famLayout;{element.FamName}, [0 0 0 0], 0];
            j=j+1;
        end

    end
    if (strcmpi(element.FamName,'GS')||strcmpi(element.FamName,'GE')||...
        strcmpi(element.FamName,'mon')||strcmpi(element.FamName,'bpm')    )
        famLayout=[famLayout;{element.FamName}, [0 0 0 0], 0];
        j=j+1;
    end
    


end

   
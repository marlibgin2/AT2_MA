function famLayout=famLayout(ACHRO)
%
% Returns a list of maget family names in eth sequence they appera in the
% input lattice cell array
%% Inputs
% ACHRO : AT2 lattice cell array
%
%% Outputs
% famLayout (NX3): cell array with two coluns. The first 
%              contains strings with names of magnet families in
%              the order they appear in the lattice. The 
%              second contains the corresponding multipolar strengths 
%              of tose elentes - the PolynomB array field of the element.
%              The third contains the bend angle of the element in degrees
% 
%% Usage examples
% fLO=famLayout(ACHRO);
%

%% History
% PFT 2024/07/17

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
                if(not(strcmp(element.FamName,famLayout{j})))
                    famLayout=[famLayout;{element.FamName},element.PolynomB,bendangle*180/pi];
                    j=j+1;
                end
            else
                famLayout=[famLayout;{element.FamName},element.PolynomB,bendangle*180/pi];
                j=j+1;
            end
        end
    end
end

   
function famLayout=famLayout(ACHRO,varargin)
%
% Returns a matlab table of element names and their properties 
% in the sequence they appear in the input lattice cell array
%% Inputs
% Mandatory argments
% ACHRO : AT2 lattice cell array
%
% Optional flags
% 'unique' : suppress output when the same famity appears in sequence
%
%% Outputs
% fL (NX7): cell array with the following columns
%     Element Name
%     Length [m]
%     Bending Angle [deg]
%     InvRho[m**-1]
%     Q[m**-2]
%     S[m**-3]
%     O[m**-4]
% 
%% Usage examples
% fL=famLayout(ACHRO);
%

%% History
% PFT 2024/07/17 : first version
% PFT 2024/07/24 : added girder markers, correctors and monitors 
%                  if available
% PFT 2024/08/03 : fixed bug update counter for correctors
% PFT 2024/08/07 : added element lengths to the output and broke up 
%                  PolynomB output into separate columns.
%                  changed output from cell array to table
%                  changed order of columns
%                  improved handling of PolynomB vector with less
%                  than 4 elements
%                  added Spos to the output table
%
%% Input argument parsing
uniquef      = any(strcmpi(varargin,'unique'));

famLayout={};
j=0;
for i=1:numel(ACHRO)
    element=ACHRO{i};
    Spos=findspos(ACHRO,i);
    if (isfield(element,'PassMethod'))
        PM=element.PassMethod;
        if (strcmp(PM,'BndMPoleSymplectic4Pass')||strcmp(PM,'StrMPoleSymplectic4Pass'))
            if (isfield(element,'BendingAngle'))
                bendangle=element.BendingAngle;
            else
                bendangle=0.0;
            end
            PolynomB=element.PolynomB;
            max=size(PolynomB,2);
            if (max<4)
                PolynomB(1,max+1:4)=0.0;
            end
            if (j>0)
                if(not(strcmp(element.FamName,famLayout{j,1}))||not(uniquef))
                    famLayout=[famLayout;Spos,{element.FamName},element.Length,...
                        bendangle*180/pi,...
                        PolynomB(2),PolynomB(3),PolynomB(4)];
                    j=j+1;
                end
            else
                famLayout=[famLayout;Spos,{element.FamName},element.Length,...
                    bendangle*180/pi,...
                    PolynomB(2),PolynomB(3),PolynomB(4)];
                j=j+1;
            end
        end
        if (strcmp(PM,'CorrectorPass'))
            famLayout=[famLayout;Spos,{element.FamName},0.0,element.Length, 0, 0, 0];
            j=j+1;
        end

    end
    if (strcmpi(element.FamName,'GS')||strcmpi(element.FamName,'GE')||...
        strcmpi(element.FamName,'mon')||strcmpi(element.FamName,'bpm')    )
        famLayout=[famLayout;Spos,{element.FamName},0.0, element.Length, 0, 0, 0];
        j=j+1;
    end
end

famLayout=cell2table(famLayout,'VariableNames', ["S[m]" "Element" "Length[m]" "Bending Angle[deg]" "Q[m**-2]" "S[m**-3]" "O[m**-4]"]);


 

   
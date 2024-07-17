function CLv = chalevel(varargin)
% Calculates the Challenge Levels of a set of magnet strengths
%   
%% Inputs
% Mandatory arguments
% Magnet Strength Limits: structure desbring magnet strength limits
%
% Optional arguments
% mode: 'X0'   - input is a vector of magnet strengths
%       'LS' - input is a lattice structure
%       default is LS
%
% X0   : (1XN) array of magnet strengths
% eqfam: (1XN) cell array of strings with 
%              the names in the MagnetStrengthLmits structure that are to
%              be compared with X0, default ={}
%
% eqsca: (1XN) array of scaling factors: useful to scale the dipole 
%         gradients before comparison, when the dipole bending angle 
%         has been changed to compensate for the introduction of 
%         reverse bends, default = ones(1,n)
% 
% ACHRO : Lattice structure from whihcv magnet stregths are to be taken
% Fams  : (1XN) cell array of strings with the names of families in the
%
% verbose :defines level of verbose output, default=0, i.e. no output
%% Outputs
% CLv : structure with fields
%   CLv.date = calculation date;
%   CLv.inputs.MagnetStrengthLimits = input Magnet Strengthstructure
%   CLv.inputs.ACHRO
%   CLv.inputs.Fams
%   CLv.inputs.X0 
%   CLv.inputs.eqfam
%   CLv.inputs.eqscal
%
% CLV.outpits.X0 : may be diferent frominput in mode LS.  
% CLv.outputs.CL   = (1xN) aray of challlenge Levels;
%
%% Usage Examples
% Fams = {'Q1_d1';'Q2_d1';'R1_d1';...
%          'D2_d1';'D3_d1';'D1_d1';...
%          'Q3_d1';'Q4_d1';'S1_d1';...
%          'S2_d1';'S3_d1';'S4_d1';...
%          'S5_d1';'O1_d1';'O2_d1';...
%          'O3_d1';'T1_d1';'T2_d1'};  
%   
%  eqfam = {'Qfend_Qdend';'Qfend_Qdend';'Qf_Qfm';...
%           'dip';      'dip';        'dipm';...
%           'Qf_Qfm';    'Qf_Qfm';     'Sdend';...
%           'Sfm';       'Sd';         'Sfi_Sfo';...
%           'Sfi_Sfo';   'Oxx_Oxy';    'Oxx_Oxy';...
%           'Oyy';       'SextTrimsQ'; 'SextTrimsQ'};
%  eqsca = [1 1 1 1 (3+4*3.49E-3*180/pi)/3 1 1 1 1 1 1 1 1 1 1 1 1 1];
%
% CLv = chalevel(MagnetStrengthLimits,'mode','X0','X0',X0,'eqfam',...
%                eqfam, 'eqsca', eqsca )';
%
% CLv = chalevel(MagnetStrengthLimits,'mode','LS','ACHRO',ACHRO_a1,'eqfam',...
%                eqfam, 'eqsca', eqsca, 'Fams', Fams )';
%

%% History
% PFT 2024/05/29: first version
% PFT 2024/06/02: documentation
% PFT 2024/06/03: further documentation
% PFT 2024/06/05: added possibility of lattice structure input
% PFT 2024/06/18: added handling of case without MagnetStrengthLimist
%                 structure (for use with the cLatt function). added verbose
%                 level parameter
% PFT 2024/07/16: added handlign of case wjere specific families have no
%                 equivaent fdamily in ManetStregtLimits, which is
%                 indicated by an empty string entry in eqfam
%
%% Input argument parsing
MagnetStrengthLimits = getargs(varargin,[]);

mode  = getoption(varargin,'mode','LS');
eqfam = getoption(varargin, 'eqfam',...
            {'Qfend_Qdend';'Qfend_Qdend';'Qf_Qfm';...
             'dip';      'dip';        'dipm';...
             'Qf_Qfm';    'Qf_Qfm';     'Sdend';...
             'Sfm';       'Sd';         'Sfi_Sfo';...
             'Sfi_Sfo';   'Oxx_Oxy';    'Oxx_Oxy';...
             'Oyy';       'SextTrimsQ'; 'SextTrimsQ'});

eqsca = getoption(varargin,'eqsca',[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
X0    = getoption(varargin,'X0', zeros(1,18));
ACHRO = getoption(varargin,'ACHRO',[]);
Fams  = getoption(varargin,'Fams', {'Q1_d1';'Q2_d1';'R1_d1';...
                           'D2_d1';'D3_d1';'D1_d1';...
                           'Q3_d1';'Q4_d1';'S1_d1';...
                           'S2_d1';'S3_d1';'S4_d1';...
                           'S5_d1';'O1_d1';'O2_d1';...
                           'O3_d1';'T1_d1';'T2_d1'} );

verboselevel = getoption(varargin,'verbose',0);
%% Extacts magnet strengths
switch mode
    case {'X0';'x0'}
        if (isempty(X0))
            fprintf('%s chalevel Error: empty X0 array for mode X0, aborting...\n',datetime);
            CLv=struct;
            return
        end
        CLv.inputs.X0=X0;

    case {'LS';'ls'}
        
        if (isempty(ACHRO)||not(iscell(ACHRO)))
            fprintf('%s chalevel Error:invalid ACHRO input for mode LS, aborting...\n',datetime);
            CLv=struct;
            return
        end
        if ( (numel(Fams)~=numel(eqfam))||(numel(Fams)~=numel(eqsca)))
            fprintf('%s chalevel Error:invalid ACHRO input for mode LS, Check that number of elements in Fams, eqfam and eqsca are the same, aborting...\n',datetime);
            CLv=struct;
            return
        end

        nfams = numel(Fams);
        X0=zeros(1,nfams);
        for i=1:nfams
            I_fam = find(atgetcells(ACHRO, 'FamName', Fams{i}));
            if (isempty(I_fam))
               fprintf('%s Error: Family %s, not found, skipping...\n',datetime, Fams{i});
               X0(i)=nan;
               continue
            end
            PassMethod = atgetfieldvalues(ACHRO, I_fam, 'PassMethod');
            K_fam =  atgetfieldvalues(ACHRO, I_fam, 'PolynomB');
            Ks = zeros(1,numel(I_fam));
            for j=1:numel(I_fam)
                switch PassMethod{1}
                    case 'BndMPoleSymplectic4Pass'
                        Ks(j)=K_fam{j}(2);
                    otherwise
                        [Kmmax,jpos] =  max(abs(K_fam{j}));
                        Ks(j) = Kmmax*sign(K_fam{j}(jpos));
                end
            end
            [Kmax, ipos] = max(abs(Ks));
             X0(i) = Kmax*sign(Ks(ipos));
        end

    otherwise
        fprintf('%s Error: Unknown mode %s , aborting...\n',datetime, mode);
end
nfams=numel(X0);

X0scal=abs(X0./eqsca);
CL=zeros(1,nfams);


%% Calculate challenge levels
if (isempty(fieldnames(MagnetStrengthLimits)))
    if (verboselevel>0)
        fprintf('%s chlevel Warning: MagnetStrengthLimits Structure not available\n',datetime);
    end
    CL=[];
else
    for i=1:nfams
        if (isfield(MagnetStrengthLimits,eqfam{i}))
            MagStr=MagnetStrengthLimits.(eqfam{i});
            Mins= MagStr.Mins;
            Maxs= MagStr.Maxs;
            CLs = MagStr.CLs; 
            ncl = numel(Mins);
            CLa = repmat(X0scal(i),1,ncl);
            CLb = (CLa>=Mins).*(CLa<Maxs);
            CLc = repmat(1000,1,ncl);
            for j=1:ncl
                if(CLb(j))
                    CLc(j)=CLs(j);
                end
            end
            CL(i)=min(CLc);
        else
            CL(i)=nan;
        end
    end
end
%% Collects output structure data
CLv.date=sprintf('%s', datetime);
CLv.inputs.MagnetStrengthLimits=MagnetStrengthLimits;
CLv.inputs.mode=mode;
CLv.inputs.eqfam=eqfam;
CLv.inputs.eqsca=eqsca;
CLv.inputs.Fams=Fams;
CLv.inputs.ACHRO=ACHRO;

CLv.outputs.X0=X0;
CLv.outputs.CL = CL;
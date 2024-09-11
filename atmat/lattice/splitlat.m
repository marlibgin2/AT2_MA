function ACHRO_SP=splitlat(ACHRO,split)
%
% generates a lattice by splitting each element 
% in the input lattice by a give factot
%
%% Inputs
% Mandatory arguments
%
% ACHRO : AT2 lattice cell array
% split  : number of splits per element
%
%% Outputs
% ACHRO_SP : splitted lattice
%
%% Usage examples
% ACHRO_SP = splitlat(ACHRO,10;

%% History
% PFT 2024/07    : first verstion
% PFT 2024/08/08 : added documentation

    splits=1/split*ones(1,split);
    ACHRO_SP={};
    for i=1:numel(ACHRO)
        if (ACHRO{i}.Length>0.0)
            element=atdivelem(ACHRO{i},splits);
            ACHRO_SP=[ACHRO_SP;element];
        else
            ACHRO_SP=[ACHRO_SP;ACHRO{i}];
        end
    end
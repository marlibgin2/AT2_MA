function [deltamax]=momentum_aperture_at(varargin)
% function [deltamin, deltamax...
%     ]=momentum_aperture_at(THERING,...
%     deltalimit,... [min max]
%     initcoord,... [x y] initial coordinate
%     delta,...
%     precdelta,...
%     deltastepsize,...
%     splits,... % number of splitting
%     split_step_divisor,...  % divide the step size at every split
%     nturns)  % Number of turns to track
%
% following the ELEGANT routine:
% Start with ? = 0, i.e., zero momentum offset.
% 2. Track a particle to see if it gets lost. If so, proceed to step 4.
% 3. Increase ? by step size ?? and return to step 2.
% 4. If no splitting steps remain, proceed to the next step. Otherwise:
% (a) Change ? to deltas ? sb??., where ?s is the largest ? for which the particle survived, and sb is the steps_back parameter.
% (b) Divide the step size by split_step_divisor to get a new step size ??.
% (c) Set?=?+??.
% (d) Decrement the ?splits remaining? counter by 1.
% (e) Continue from step 2.
% 5. Stop. The momentum aperture is ?s
%
% ex: [deltamax]=momentum_aperture_at(THERING,0.1,[10^-6 10^-6],0,0,0.01,3,10,100)
% ex: [deltamin]=momentum_aperture_at(THERING,-0.1,[10^-6 10^-6],0,0,-0.01,3,10,100)

%disp([delta splits])

% <MSj, MAX IV 2024-02-22> Due to the recursive function calls, significant
% slowdowns can be expected if the amount of inputs are not minimized.
% Rewrite to establish persistent data structures in the first function
% calls.

persistent THERING deltalimit initcoord split_step_divisor nturns

%% Input handling

if ~iscaller('momentum_aperture_at')
    % On first call, set all the static information
    THERING = varargin{1};
    deltalimit = varargin{2};
    initcoord = varargin{3};
    delta = varargin{4};
    precdelta = varargin{5};
    deltastepsize = varargin{6};
    splits = varargin{7};
    split_step_divisor = varargin{8};
    nturns = varargin{9};

    % Track one turn to establish the data structures in atpass memory
    ringpass(THERING,[initcoord(1) 0 initcoord(2) 0 delta 0]',1);
else
    % For recursive calls, only changing information should have been
    % passed...
    delta = varargin{1};
    precdelta = varargin{2};
    deltastepsize = varargin{3};
    splits = varargin{4};
end

%% Recursive search
if ( delta>=0 && delta<deltalimit) ||  ( delta<=0 && delta>deltalimit)
        
    if splits>-1
                
        % Track for this delta
        %  Faster if tracking several particles per call? Not clear,
        %  depends on overhead and parallellization in atpass.
        [~, LOSS] =ringpass(THERING,[initcoord(1) 0 initcoord(2) 0 delta 0]',nturns,'KeepLattice');
        
        if LOSS~=1 % if NOT LOST go to next step
            [deltamax...
                ]=momentum_aperture_at(delta+deltastepsize,... % delta center
                delta,...
                deltastepsize,...
                splits);
            
        else % if LOST reduce stepsize
            [deltamax...
                ]=momentum_aperture_at(precdelta+deltastepsize/split_step_divisor,... % go back to previous delta center and increase of smaller step
                precdelta,...
                deltastepsize/split_step_divisor,...
                splits-1);
        end
    else
        
        % no splitting steps remain
        deltamax=delta-deltastepsize;
        
    end
    
else
    % limit reached
    deltamax=delta;
end


return;




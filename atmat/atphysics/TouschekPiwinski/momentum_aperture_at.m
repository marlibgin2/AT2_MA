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

% <MSj, MAX IV 2024-02-22> A number of alterations done:
%   a) Due to the recursive function calls, significant slowdowns can be
%      expected if the amount of inputs are not minimized. I.e. passing a
%      large lattice an amount of times corresponding to the 'splits'
%      parameter can be costly. To avoid this the function was slightly
%      rewritten to establish persistent data structures in the first
%      function calls. 
%   b) As the lattice is not altered the ringpass flag 'KeepLattice' was
%      used. NB! This has been verified to work even in parfor loops; each
%      worker has it's own set of persistent variables in MATLAB. 
%   c) Input handling altered to provide reasonable default inputs if any
%      parameter is missing or empty. In particular:
%        deltalimit = 0.1
%        splits = 10
%        split_step_divisor = 2
%      This results in a final precision around 0.1/(2^10), i.e. 1e-4.
%      Furthermore:
%        initcoord = findorbit6(THERING,1)
%        nturns = <full synchrotron period, as calculated by atx>
%      Neither of the above are particularly time consuming compared to the
%      actual tracking calculations. The starting orbit does make a
%      slight difference to the result based on empirical tests, although
%      this may not be that surprising.

persistent THERING deltalimit initcoord split_step_divisor nturns

%% Input handling

if ~isrecursivecall
    % On first call, set all the static information
    % If the information haven't been given for some parameters, react
    % accordingly with either an error or default settings

    if nargin < 1 || isempty(varargin{1})
        error('momentum_aperture_at:No lattice given!');
    else
        THERING = varargin{1};
    end

    if nargin < 2 || isempty(varargin{2})
        % Assume upper limit of 10%
        warning('momentum_aperture_at:No upper search boundary defined, using +-10%.');
        deltalimit = 0.1;
    else
        deltalimit = varargin{2};
    end

    if nargin < 3 || isempty(varargin{3})
        % Calculate closed orbit and use as input.
        if check_6d(THERING)
            warning('momentum_aperture_at:No initial orbit given, investigating using momentum deviations from 6D closed orbit.');
            initcoord = findorbit6(THERING,1);
        else
            warning('momentum_aperture_at:No initial orbit given, investigating using momentum deviations from 4D closed orbit.');
            initcoord = [findorbit4(THERING,1);0;0];
        end
    else
        initcoord = varargin{3};
    end
    % Depending on the input dimensions, produce the full 6D vector
    if numel(initcoord) == 2
        initcoord = [initcoord(1); 0; initcoord(2); 0; 0; 0];
    elseif numel(initcoord) == 4
        initcoord = [initcoord; 0; 0];
    end


    if nargin < 4 || isempty(varargin{4})
        delta = deltalimit / 2;
    else
        delta = varargin{4};
    end

    if nargin < 5 || isempty(varargin{5})
        precdelta = 0;
    else
        precdelta = varargin{5};
    end

    if nargin < 6 || isempty(varargin{6})
        deltastepsize = delta / 2;
    else
        deltastepsize = varargin{6};
    end

    if nargin < 7 || isempty(varargin{7})
        splits = 10;    % NB! This together with a split_step_divisor of 2 corresponds to a precision of deltalimit / 1024. So 4th digit.
    else
        splits = varargin{7};
    end

    if nargin < 8 || isempty(varargin{8})
        split_step_divisor = 2;
    else
        % Question here is whether anything else than 2 is a sane choice?
        split_step_divisor = varargin{8};
    end

    if nargin < 9 || isempty(varargin{9})
        % If the amount of turns are not specified, use slightly more than
        % one synchrotron period (1.2). 
        % NB! For some reason there are slight differences between 2 and 1
        % synchrotron period even though the particle should have gone
        % through the momentum extremes.
        [~, params] = atx(THERING);
        nturns = round(1.2 * (1/params.fs) * PhysConstant.speed_of_light_in_vacuum.value / params.ll);
    else
        nturns = varargin{9};
    end
    

    % Track one turn to establish the data structures in atpass memory
    ringpass(THERING,initcoord,1);
else
    % For recursive calls, we have more control of the input. Only changing
    % information should be passed...
    delta = varargin{1};
    precdelta = varargin{2};
    deltastepsize = varargin{3};
    splits = varargin{4};
end

%% Notes for possible improvement
% Without parallellizing the calculations it is worth noting that, due to
% some overhead, it's possible to track three particles in less time it
% takes to track one particle twice. Tracking four particles takes roughly
% as much time as tracking one particle twice. Hence a speed-up is likely
% possible if one swaps to track three particles, distributed so that the
% search area is split into four parts. Estimates indicate this is only on
% the 10% level though, so nothing massive.
%

%% Recursive search
if ( delta>=0 && delta<=deltalimit) || ( delta<=0 && delta>=deltalimit)
        
    if splits>-1
                
        % Track for this delta
        %  Note that a slight transverse orbit shift is added based on
        %  discussions on the elegant forum.
        [~, LOSS] = ringpass(THERING,initcoord + [10e-6 0 10e-6 0 delta 0]',nturns,'KeepLattice');
        
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

end


function recursive = isrecursivecall
stack=dbstack;
if numel(stack) < 3
    recursive = false;
else
    recursive = strcmp(stack(2).name,stack(3).name);
end
end

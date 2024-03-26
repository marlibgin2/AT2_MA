function [map_l,map_h, varargout]=calculateLMA(RING,varargin)
%CALCULATELMA computes local momentum aperture for the first lattice period
% [map_l,map_h]=calculateLMA(RING,points)
% [map_l,map_h]=calculateLMA(...,PARAMETER,VALUE,...)
%
% INPUT
% 1. RING       {cell array of structs, required} AT2 lattice for which the
%               LMA should be calculated
% 2. points     {column vector} index elements for which the local momentum
%               aperture should be calculated
% 3. PARAMETER, VALUE   These pairs are passed onto momentum_aperture_at
%                       and corresponds to the standard inputs to that
%                       function. 
%                       
% NOTES
% 1. In order to determine the first lattice period/achromat, the function
%    first searches for a lattice marker element with 'FamName' equal to
%    'AchrEnd'. If this is not found it calls estimateperiodicity to get a
%    best guess.
% 2. The routine makes use of parfor in order to speed things up.
% 3. Please ensure the input lattice has been configured for either 6D or
%    4D, depending on what you wish to investigate.
%
% See also momentum_aperture_at, estimateperiodicity

% QUESTION MARKS 
% a) How granular does the lattice have to be? First rough run, then split
% elements if necessary?
%
% all Ring momentum aperture
% points=1:10:length(THERING);

points = [];
if nargin > 1 && isnumeric(varargin{1})
        points = varargin{1};
        varargin = varargin(2:end);
end

% Attempt to restrict to a single achromat
if isempty(points)
    % First check whether the first achromat has been explicitly marked...
    I = findcells(RING,'FamName','AchrEnd');
    if isempty(I)
        % If not, estimate the periodicity and grab the first achromat based on
        % that. NB! If there are superperiods the below will flag the first
        % superperiod to be calculated.
        [~, I] = estimateperiodicity(RING);
        I = I{1}(end);
    end
    J = cellfun(@(x) x.Length > 0,RING(1:I(1)));
    I = 1:I(1);
    points = I(J)';
end

map_h=zeros(length(points),1);
map_l=zeros(length(points),1);

%% Input argument parsing

% Ensure the lattice is a column vector
RING = RING(:);

% For the amount of turns to track, use an amount equivalent to 1.2
% synchrotron periods.
% NB! This is the default in momentum_aperture_at, but it's a common
% parameter and so can be calculated here for efficiency. This also ensures
% a sufficient amount of turns are used in case the RF acceptance is
% modified.
[~, params] = atx(RING);
nturns = round(1.2 * (1/params.fs) * PhysConstant.speed_of_light_in_vacuum.value / params.ll);

% Check for input flags to override default momentum_aperture_at settings.
% Empty values uses momentum_aperture_at default values.
deltalimit    = getoption(varargin,'deltalimit',0.1);
delta         = getoption(varargin,'delta',0.0);
precdelta     = getoption(varargin,'precdelta',0.0);
deltastepsize = getoption(varargin,'deltastepsize',0.05);
splits        = getoption(varargin,'splits',10);
split_step_divisor = getoption(varargin,'split_step_divisor',2);
nturns        = getoption(varargin,'nturns',nturns);
initcoord     = getoption(varargin,'initcoord',[]);


%% Calculates Momentum Aperture
% Could get a bit of performance increase by parallelizing the positive and
% negative function calls, assuming sufficient cores are available (2*nbr
% of elements). 
parfor i=1:length(points)
   % disp([i length(points) i/length(points)*100])
    %cycle ring
     RING_cycled = circshift(RING,1-points(i),1);
%      RING_cycled=[RING(points(i):end); RING(1:points(i)-1)]';
        try
            map_h(i)=momentum_aperture_at(RING_cycled,abs(deltalimit),initcoord,...
                delta,precdelta,abs(deltastepsize),splits,split_step_divisor,nturns);
        catch
            map_h(i)=0;
        end
        
        try
            map_l(i)=momentum_aperture_at(RING_cycled,-abs(deltalimit),initcoord,...
                delta,precdelta,-abs(deltastepsize),splits,split_step_divisor,nturns);
        catch
            map_l(i)=0;
        end
       
end

if nargout > 2
    varargout{1} = points;
end

end




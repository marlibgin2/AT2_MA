function values = atgetfieldvalues(ring,varargin)
%ATGETFIELDVALUES retrieves the field values AT cell array of elements
%
% VALUES = ATGETFIELDVALUES(RING,'field') extracts the values of
% the field 'field' in all the elements of RING
%
% VALUES = ATGETFIELDVALUES(RING,INDEX,'field') extracts the values of
% the field 'field' in the elements of RING selected by INDEX
%
% if RING{I}.FIELD is a numeric scalar
%    VALUES is a length(INDEX) x 1 array
% otherwise
%    VALUES is a length(INDEX) x 1 cell array
%
%VALUES = ATGETFIELDVALUES(...,'Default',default_value) Uses default_values
%   if the required field is not existing
%
% More generally ATGETFIELDVALUES(RING,INDEX,subs1,subs2,...) will call
%  GETFIELD(RING{I},subs1,subs2,...) for I in INDEX
%
% Examples:
%
% V=ATGETFIELDVALUES(RING,1:10,'PolynomB') is a 10x1  cell array
% such that V{I}=RING{I}.PolynomB for I=1:10
%
% V=ATGETFIELDVALUES(RING(1:10),'PolynomB',{1,2}) is a 10x1 array
% such that V(I)=RING{I},PolynomB(1,2)
%
%
% See also ATSETFIELDVALUES ATGETCELLS GETCELLSTRUCT FINDCELLS

[default_val,vargs]=getoption(varargin,'Default',NaN);
def=(length(vargs)==length(varargin));  % No default value specified

if islogical(vargs{1}) || isnumeric(vargs{1})
    index = vargs{1}; vargs = vargs(2:end);
end

% if numel(vargs) > 1
%     pos = vargs{2:end};
% else
%     pos = {':'};
% end

        if exist('index','var')
            [values,isnumscal,isok]=cellfun(@scan,ring(index),'UniformOutput',false);
        else
            [values,isnumscal,isok]=cellfun(@scan,ring,'UniformOutput',false);
        end
        isok=cell2mat(isok);
        isnumscal=cell2mat(isnumscal);
        if all(isnumscal)
            values=cell2mat(values);
        elseif all(~isnumscal(isok)) && def
            values(~isok)={[]};
        end

        % NB! Any checks, filtering, etc. that can be moved out from this
        % function the better as it will be called a LOT!
        function [val,isnumscal,isok]=scan(el)
            if isfield(el,vargs{1})
%                     val = el.(vargs{1})(vargs{2:end});  % NB! Calling with dynamic field names is faster but less flexible.
                val=getfield(el,vargs{:});   
                isok=true;
            else
                val=default_val;
                isok=false;
            end

% NB! try/catch statements are expensive and are to be avoided, especially
% in low-level functions likely to be called thousands of times!
%             try
%                 val=getfield(el,varargin{:});
%                 isok=true;
%             catch
%                 val=default_val;
%                 isok=false;
%             end
            isnumscal=isnumeric(val) && isscalar(val);
        end
   

end

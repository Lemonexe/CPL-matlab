%s = EXTENDSTRUCT(varargin)
%       extendStruct creates struct 's' from parameters just like struct() does
%       extendStruct can, however, accept a struct as first parameter and extend it with the following parameters
%       properties of the initial struct may be overwritten with input parameters
%
%       s0 = extendStruct('a', 1, 'b', 2)
%       s1 = extendStruct(s0, 'b', 3, 'c', 3, 'd', 4)

function s = extendStruct(varargin)
    if nargin > 0 && isstruct(varargin{1})
        s = varargin{1};
        if mod(numel(varargin), 2) == 0; error('When first argument is a struct, the other field and value input arguments must come in pairs.'); end
        for i = 2:2:(numel(varargin)-1)
            s.(varargin{i}) = varargin{i+1};
        end
    else
        if mod(numel(varargin), 2) == 1; error('Field and value input arguments must come in pairs.'); end
        s = struct();
        for i = 1:2:numel(varargin);
            s.(varargin{i}) = varargin{i+1};
        end
    end
end

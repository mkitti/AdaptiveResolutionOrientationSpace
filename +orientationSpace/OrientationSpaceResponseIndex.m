classdef OrientationSpaceResponseIndex
    %orietationSpace.OrientationSpaceResponseIndex allow for indexing into a
    %OrientationResponse object
    
    properties
        response
    end
    
    methods
        function obj =  OrientationSpaceResponseIndex(orientationSpaceResponse)
            obj.response = orientationSpaceResponse;
        end
        function varargout = subsref(obj,S)
            switch(S(1).type)
%                 case '.'
                case '()'
                    varargout{1} = obj.response.getResponseAtPoint(S(1).subs{:});
                case '{}'
                    varargout{1} =  obj.response.getResponseAtOrientation(S(1).subs{:});
                otherwise
                    [varargout{1:nargout}] = builtin('subsref',obj,S);
            end
        end
    end
    
end


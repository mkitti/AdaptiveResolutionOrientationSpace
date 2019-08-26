classdef OrientationSpaceNMS < handle
    %OrientationSpaceNMS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        response
    end
    properties (Dependent = true)
        nms;
    end
    
    methods
        function obj = OrientationSpaceNMS(response)
            obj.response = response;
        end
        function h = imshow(obj,varargin)
            h = imshow(obj.nms,varargin{:});
        end
        function h = imshowpair(A,B)
            if(nargin > 1)
                if(isa(B,'OrientationSpaceNMS'))
                    B = real(B.nms);
                end
            else
                B = imag(A.nms);
            end
            if(isa(A,'OrientationSpaceNMS'))
                A = real(A.nms);
            end
            h = imshowpair(A,B);
        end
        function nms = get.nms(obj)
            nms = obj.response.nms;
        end
    end
    
end


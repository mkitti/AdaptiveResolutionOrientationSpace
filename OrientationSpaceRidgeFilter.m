classdef OrientationSpaceRidgeFilter < OrientationSpaceFilter
    %OrientationSpaceRidgeFilter OrientationSpace filter that does line
    %filtering only
    
    properties
    end
    
    methods
        function obj = OrientationSpaceRidgeFilter(f_c,b_f,K,varargin)
            obj@OrientationSpaceFilter(f_c,b_f,K,varargin{:});
        end
        function R = getResponse(obj,I)
            If = fft2(I);
            angularResponse = obj.applyRidgeFilter(If);
            ns = [0 cumsum([obj.n])];
            R(numel(obj)) = OrientationSpaceResponse;
            for o=1:numel(obj)
                R(o) = OrientationSpaceResponse(obj(o),angularResponse(:,:,ns(o)+1:ns(o+1)));
            end
            R = reshape(R,size(obj));
        end
    end
    
end


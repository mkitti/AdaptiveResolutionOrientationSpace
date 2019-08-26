%nms = nonMaximumSuppression(res, th, sup)
%
% Inputs:   res : response
%            th : orientation
%           sup : value to set suppressed pixels to
%
% Uses the grid conventions of steerableDetector()

% Francois Aguet

function res = nonMaximumSuppression(res, th, sup)

if(nargin < 3)
    sup = 0;
end

[ny,nx] = size(res);

res = padarrayXT(res, [1 1], 'symmetric');

[x,y] = meshgrid(1:nx,1:ny);

% +1 interp
A1 = interp2(res, x+1+cos(th), y+1+sin(th),'linear',0);

% -1 interp
A2 = interp2(res, x+1-cos(th), y+1-sin(th),'linear',0);

res = res(2:end-1,2:end-1);

res(res<A1 | res<A2) = sup;

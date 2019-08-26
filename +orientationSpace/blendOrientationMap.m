function [ blended_map ] = blendOrientationMap( theta, res, cm )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 2 || isempty(res))
    res = ones(size(theta));
end
if(nargin < 3)
    cm = hsv(256);
end

theta_not_nan = ~isnan(theta);
theta_mapped = NaN([numel(theta) 3]);
theta_mapped(theta_not_nan,:) = cm(floor(theta(theta_not_nan)/pi*(size(cm,1)-1))+1,:);
sz_theta = size(theta);
sz_theta(ndims(theta)+1:3) = 1;
theta_mapped = reshape(theta_mapped,[sz_theta 3]);
res(res < 0) = NaN;
blended_map = squeeze((nansum(theta_mapped.*res,3)./nansum(res,3)));

if(nargout == 0)
    figure;
    imshow(blended_map);
end


% figure;
% imshow(squeeze((nansum(theta_mapped.*res,3)./nansum(res,3))),[]);

end


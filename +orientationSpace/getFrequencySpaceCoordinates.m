function [ polar, cartesian ] = getFrequencySpaceCoordinates( N )
%GETFREQUENCYSPACECOORDINATES orientationSpace.getFrequencySpaceCoordinates
% obtain polar and cartesian coordinates
%
% Caches the last coordinates for the last N
%
% INPUT
% N - size of space, if scalar a size of [N N] is extrapolated
%
% OUTPUT
% polar - polar coordinates, shifted
%    .f - frequency magnitude, from 0 to sqrt(2)/2
%    .theta - frequency angle (-pi,pi]
% cartesian - cartesian coordinates, shifted
%    .X - X coordinate from -floor(N/2) to N-floor(N/2)
%    .Y - Y coordinate from -floor(N/2) to N-floor(N/2)
% 
% See also ifftshift, fftshift

persistent lastN lastCartesian lastPolar


if(nargin < 1)
    N = 1024;
end

if(isscalar(N))
    N = [N N];
end

% if(~isempty(lastN) && all(N == lastN))
if(isequal(lastN,N))
    cartesian = lastCartesian;
    polar = lastPolar;
else   
    lastN = N;

    % [X, Y] = meshgrid(-N:(N-1));
    [cartesian.X, cartesian.Y] = meshgrid( (0:N(2)-1) - floor(N(2)/2), (0:N(1)-1) - floor(N(1)/2) );

    cartesian.X = ifftshift(cartesian.X);
    cartesian.Y = ifftshift(cartesian.Y);

    [polar.theta, polar.f] = cart2pol(cartesian.X / floor(N(2)/2) / 2,cartesian.Y / floor(N(1)/2) / 2);
    
    lastCartesian = cartesian;
    lastPolar = polar;

    % polar.theta = ifftshift(polar.theta);
    % polar.f = ifftshift(polar.f);
end


end


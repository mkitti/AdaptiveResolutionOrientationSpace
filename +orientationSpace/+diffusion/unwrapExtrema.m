function [ extrema ] = unwrapExtrema( extrema,events, period)
%unwrapExtrema Unwrap aligned extrema around the periodic boundary to avoid
%discontinuities

if(nargin < 3)
    period = 2*pi;
end

change = extrema(:,events+1) - extrema(:,events);
[r,c] = find(abs(change) > period/2);

for e = 1:length(r)
    extrema(r(e),events(c(e))+1:end) = extrema(r(e),events(c(e))+1:end) -sign(change(r(e),c(e)))*period;
end


end


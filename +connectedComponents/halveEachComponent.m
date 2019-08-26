function [ cc ] = halveEachComponent( cc, minLength )
%halveEachComponent Split each component into two components at the
%midpoint
%
% Use orderPixelIdxList first

if(nargin < 2)
    minLength = 5;
end

nElements = cellfun('prodofsize',cc.PixelIdxList);
splitMe = nElements >= minLength;
halved = arrayfun(@(x,n) {x{1}(1:n) x{1}(n+1:end)},cc.PixelIdxList(splitMe),ceil(nElements(splitMe)/2),'UniformOutput',false);
if(isempty(halved))
    return;
else
    halved = [halved{:} cc.PixelIdxList{~splitMe}];
end

cc.PixelIdxList = halved;
cc.NumObjects = length(halved);


end


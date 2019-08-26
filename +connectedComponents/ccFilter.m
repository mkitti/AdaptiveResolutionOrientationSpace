function [ cc, oldccidx ] = ccFilter( cc, filter )
%ccFilter Filter a connected components structure
%
% cc is a connected components structure from bwconncomp
% filter is a vector representing a subindex into PixelIdxList
%        used to select the new connected components. 
%        This may be a logical index or one selecting specific
%        components of the connected components structure
%
% output:
%   newcc is a new conneccted component struccture with filtere applied
%   oldccidx is a vector containing numeric indices of the original
%     connected component struccture. Useful if logical indexing was used.

% newcc = cc;
cc.PixelIdxList = cc.PixelIdxList(filter);
cc.NumObjects = length(cc.PixelIdxList);
if(nargout > 1)
    idx = 1:cc.NumObjects;
    oldccidx = idx(filter);
end

end


function [ bw, count, cumBranchPoints ] = bwRemoveBranchPoints( bw, n )
%BWREMOVEBRANCHPOINTS Remove all branch points
%
% INPUT
% bw - logical image, any non-logical image will be converted using logical
% n  - (optional) Number of times to remove branchpoints using bwmorph
%      default: Inf (until no more branchpoints are found)
%
% OUTPUT
% bw
% count
% branchPoints

if(nargin < 2)
    n = Inf;
end

done = false;
count = 0;
if(nargout > 2)
    cumBranchPoints = false(size(bw));
end

while(~done)
    branchPoints = bwmorph(bw,'branchpoints');
    bw = bw & ~branchPoints;
    done = ~any(branchPoints(:));
    count = count + 1;
    if(nargout > 2)
        cumBranchPoints = cumBranchPoints | branchPoints;
    end
    if(count > n)
        break;
    end
end


end


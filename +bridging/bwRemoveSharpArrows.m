function [ bw, sharpArrows ] = bwRemoveSharpArrows( bw )
%bwRemoveSharpArrows Remove 8-connected sharp arrows
%
% INPUT
% bw - binary image
%
% OUTPUT
% bw - binary image without sharp arrows
% sharpArrows - location of sharp arrows
%
% Sharp arrows are
% [ 1 0 0
%   0 1 0
%   1 0 0]
% [ 1 0 1
%   0 1 0
%   0 0 0]
% [ 0 0 0
%   0 1 0
%   1 0 1]
% [ 0 0 1
%   0 1 0
%   0 0 1]

dbw = double(logical(bw));
neighborCount = imfilter(dbw,[1 1 1; 1 0 1; 1 1 1]);
forwardSlashCount = imfilter(dbw,eye(3));
backwardSlashCount = imfilter(dbw,rot90(eye(3)));
xCount = imfilter(dbw,[1 0 1; 0 0 0; 1 0 1]);

sharpArrows = neighborCount == 2 & forwardSlashCount ~= 3 & backwardSlashCount ~= 3 & xCount == 2;
bw(sharpArrows) = 0;

end


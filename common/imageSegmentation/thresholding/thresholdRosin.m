function level = thresholdRosin(imageIn,varargin)
% Select a thresholding level using Rosin's method
%
% level = thresholdRosin(imageIn)
% level = thresholdRosin(imageIn,showPlots)
% 
% This function selects a threshold for the input fluorescence image by
% analyzing the image's intensity distribution. This requires good signal-
% to-noise, and a significant amount of background in the image. The
% threshold is selected using Rosin's method
% 
% Input:
% 
%   imageIn - The N-Dimensional image to be thresholded.
% 
% 
%   showPlots - If true, a plot of the histogram and an overlay of the mask
%   on the image will be shown. The overlay plot only works if the image is
%   2D.
% 
% 
% Output:
% 
% 
%   level - The intensity value selected for thresholding.
%
% Revamped from RosinSeg
%
% Sebastien Besson, 5/2011

ip=inputParser;
ip.addRequired('imageIn',@isnumeric);
ip.addOptional('showPlots',0,@isnumeric)
ip.parse(imageIn,varargin{:});
showPlots=ip.Results.showPlots;

%Convert to double if necessary
imageIn = double(imageIn);

%find nonzero values (due to masking)
nzInd = find(imageIn);

%get minumum and maximum pixel values in image
minSignal = min(imageIn(nzInd));
maxSignal = max(imageIn(nzInd));

%normalize nonzero value between 0 and 1
imageInNorm = zeros(size(imageIn));
imageInNorm(nzInd) = (imageIn(nzInd)- minSignal) / (maxSignal - minSignal);

[~,level] =  cutFirstHistMode(imageInNorm,0);

level = level*(maxSignal - minSignal)+minSignal;

if showPlots 
    imageMask = imageIn >= level;
    figure;
    imagesc(imageIn);
    hold on
    contour(imageMask,'w')
    colormap hot
end

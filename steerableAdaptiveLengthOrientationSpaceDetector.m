function [ response, theta, nms, angularResponse, other ] = steerableAdaptiveLengthOrientationSpaceDetector( I, varargin )
%steerableAdaptiveLengthOrientationSpaceDetector Perform adaptive length orientation
%space segmentation
%
% INPUT
% I - image
%     Type: 2D numeric matrix, non-empty
% order - (optional), K parameter for OrientationSpaceFilter
%     Type: Numeric, scalar
%     Default: 8
% sigma - (optional), scale parameter setting the radial bandpass in pixels
%         central frequency of the bandpass fliter will be 1/2/pi/sigma
%     Type: Numeric, scalar
%     Default: 2 (pixels)
%
% PARAMETERS
% filter - OrientationSpaceFilter to use, overrides order and sigma
% response - OrientationSpaceResponse to use, overrides order and sigma
% adaptLengthInRegime - Search current maxima regime to find smallest slopes
%     Type: logical
%     Default: true
% meanThresholdMethod - function to determine threshold of mean response
%     Type: char, function_handle
%     Default: @thresholdOtsu
% mask - binary mask the same size as I
%     Type: logical
%     Default: []
% nlmsMask - binary mask the same size as I to apply to NLMS
%     Type: logical
%     Default: []
% nlmsThreshold - threshold to apply to NLMS
%     Type: numeric, 2D
%     Default: []
% useParallelPool - true if parallel pool should be used
%     Type: logical
%     Default: true
% maskDilationDiskRadius - radius in pixels to dilate the mask 
%     Type: numeric
%     Default: 3
% maskFillHoles - fill holes in the mask
%     Type: logical
%     Default: false
% diagnosticMode - true if diagnostic figures should be shown
%     Type: logical, scalar
%     Default: false
% K_sampling_delta - interval to sample K when using adaptLengthInRegime
%     Type: numeric, scalar
%     Default: 0.1
% responseOrder - order from which to retrieve response (not maxima)
%     Type: numeric, scalar
%     Default: 3
%
% OUTPUT
% response - response values (at responseOrder) of maxima
%            size(image) x M
% theta    - maxima
% nms      - non-maximum suppression like image
%            same size as image
% angularResponse
%          - response values at regular intervals
% other    - struct containing various variables:
% .nlms_highest - NLMS using highest regime maxima
% .maxima_highest - Maxima at highest regime
% .nlms_highest_mip - Maximum Intensity Projection of .nlms_highest
% .maximum_single_angle - Angle at lowest single-orientation regime
% .nlms_single - NLMS at maximum_single_angle
% .nlms_single_binary - binarized .nlms_single by meanResponse
% .attenuatedMeanResponse - meanResponse attenuated by neighborhood
% .meanResponse - mean orientationSpace response
% .nlmsMask - mask used to process NLMS
% .params - input parameters
% .combinedResR - currently the same as response
%
% Mark Kittisopikul, Ph.D.
% Goldman Lab
% Northwestern University
% 2018
    
ip = inputParser;
ip.addRequired('I',@(x) validateattributes(x,{'numeric'},{'2d','nonempty'}));
ip.addOptional('order',8,@(x) validateattributes(x,{'numeric'},{'scalar'}));
ip.addOptional('sigma',2,@(x) validateattributes(x,{'numeric'},{'scalar'}));
ip.addParamValue('filter',[],@(x) validateattributes(x,{'OrientationSpaceFilter'},{'scalar'}));
ip.addParamValue('response',[],@(x) validateattributes(x,{'OrientationSpaceResponse'},{'scalar'}));
ip.addParamValue('adaptLengthInRegime',true,@islogical);
ip.addParamValue('meanThresholdMethod','otsu',@(x) validateattributes(x,{'char','function_handle'},{})); % @(x) validatestring(x,{'otsu','rosin'}));
ip.addParamValue('meanThreshold',[],@(x) validateattributes(x,{'numeric'},{'scalar'}));
ip.addParamValue('mask',[],@islogical);
ip.addParamValue('nlmsMask',[],@islogical);
ip.addParamValue('nlmsThreshold',[],@(x) validateattributes(x,{'numeric'},{'2d'}));
ip.addParamValue('useParallelPool',true,@islogical);
ip.addParamValue('maskDilationDiskRadius',3,@(x) validateattributes(x,{'numeric'},{'scalar'}));
ip.addParamValue('maskFillHoles',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
ip.addParamValue('diagnosticMode',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
ip.addParamValue('K_sampling_delta',0.1,@(x) validateattributes(x,{'numeric'},{'scalar'}));
ip.addParamValue('responseOrder',3,@(x) validateattributes(x,{'numeric'},{'scalar'}));

ip.parse(I,varargin{:});

if( isempty(ip.Results.response) )
    if( isempty(ip.Results.filter) )
        F = OrientationSpaceFilter.constructByRadialOrder(1/2/pi./ip.Results.sigma,1,ip.Results.order,'none');
    else
        F = ip.Results.filter;
    end

    R = F*I;
else
    R = ip.Results.response;
    F = ip.Results.filter;
end

if(ip.Results.useParallelPool)
    pool = gcp;
end

% Obtain Fourier transform
a_hat = fft(real(R.a),[],3);

%% Evaluate and threshold mean first to quickly reduce the problem
meanResponse = a_hat(:,:,1)./size(a_hat,3);
if(ip.Results.diagnosticMode)
    oldPrefImshowBorder = iptgetpref('ImshowBorder');
    iptsetpref('ImshowBorder','loose');
    other.diagFig(1) = figure('Name','meanResponse');
    imshow(meanResponse,[]);
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none')
end

%% Thresholding
if(isempty(ip.Results.meanThreshold))
    % If method is char, convert to a function
    if(ischar(ip.Results.meanThresholdMethod))
        switch(ip.Results.meanThresholdMethod)
            case 'otsu'
                meanThresholdMethod = @thresholdOtsu;
            case 'rosin'
                meanThresholdMethod = @thresholdRosin;
            otherwise
                meanThresholdMethod = str2func(ip.Results.meanThresholdMethod);
        end
    else
        % Assume this is a function
        meanThresholdMethod = ip.Results.meanThresholdMethod;
    end

    % We could mask the thresholding function, but ...
%     if(~isempty(ip.Results.mask))
%         meanThreshold = meanThresholdMethod(meanResponse(ip.Results.mask));
%     else
%         meanThreshold = meanThresholdMethod(meanResponse);
%     end

    % It is more flexible to leave it up to the user
    meanThreshold = meanThresholdMethod(meanResponse);
else
    % User can also set a fixed threshold directly
    % This can a single threshold or the size of the image
    meanThreshold = ip.Results.meanThreshold;
end

if(ip.Results.diagnosticMode)
    other.diagFig(end+1) = figure('Name','meanResponse Histogram');
    histogram(meanResponse(:));
    xlabel('meanResponse');
    ylabel('count');
    line([1 1].*meanThreshold,ylim,'Color','r','LineWidth',2);
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none')
end

%% Masking

if(isempty(ip.Results.nlmsMask))
    meanResponseMask = meanResponse > meanThreshold;
    if(ip.Results.diagnosticMode)
        other.diagFig(end+1) = figure('Name','meanResponse > meanThreshold');
        imshow(meanResponseMask,[]);
        title(other.diagFig(end).Children,other.diagFig(end).Name, ...
            'Interpreter','none')
    end
    if(~isempty(ip.Results.mask))
        meanResponseMask = meanResponseMask & ip.Results.mask;
    end
    if(ip.Results.diagnosticMode)
        other.diagFig(end+1) = figure('Name','meanResponseMask');
        imshow(meanResponseMask,[]);
        title(other.diagFig(end).Children,other.diagFig(end).Name, ...
            'Interpreter','none')
    end

    % For NLMS (nlmsMask)
    meanResponseMaskDilated = imdilate(meanResponseMask,strel('disk',ip.Results.maskDilationDiskRadius));
    if(ip.Results.maskFillHoles)
        meanResponseMaskDilated = imfill(meanResponseMaskDilated,'holes');
    end
    if(~isempty(ip.Results.mask))
        meanResponseMaskDilated = meanResponseMaskDilated & ip.Results.mask;
    end
    if(ip.Results.diagnosticMode)
        other.diagFig(end+1) = figure('Name','meanResponseMaskDilated');
        imshow(meanResponseMaskDilated,[]);
        title(other.diagFig(end).Children,other.diagFig(end).Name, ...
            'Interpreter','none')
        diag_rp = regionprops(meanResponseMaskDilated,'BoundingBox','Area');
        [~,diag_rp_max_idx] = max([diag_rp.Area]);
        diag_rp = diag_rp(diag_rp_max_idx);
        rectangle('Position',diag_rp.BoundingBox,'EdgeColor','r');
    end

    nlmsMask = meanResponseMaskDilated;
else
    % User defined nlmsMask
    nlmsMask = ip.Results.nlmsMask;
    if(ip.Results.diagnosticMode)
        diag_rp = regionprops(nlmsMask,'BoundingBox','Area');
    end
end

%% Setup orientation analysis problem
nanTemplate = NaN(size(nlmsMask));

a_hat = shiftdim(a_hat,2);
a_hat = a_hat(:,nlmsMask);

%% Evaluate single orientation, fast easy case

R_res = R.getResponseAtOrderFT(ip.Results.responseOrder,2);
maximum_single_angle = nanTemplate;
maximum_single_angle(nlmsMask) = wraparoundN(-angle(a_hat(2,:))/2,0,pi);
nlms_single = nonLocalMaximaSuppressionPrecise(real(R_res.a),maximum_single_angle,[],[],nlmsMask);
nlms_single_binary = nlms_single > meanResponse;

if(ip.Results.diagnosticMode)
    other.diagFig(end+1) = figure('Name','maximum_single_angle');
    
    cm = colorcet('CBC2','N',256);
    
    
    maximum_single_angle_map = orientationSpace.blendOrientationMap(maximum_single_angle,R_res.interpft1(maximum_single_angle),cm);
    imshow(maximum_single_angle_map,[]);
    caxis([-1/256 pi]);
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none')
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
    
    other.diagFig(end+1) = figure('Name','nlms_single');
    imshow(nlms_single,[]);
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none')
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
    
   other.diagFig(end+1) = figure('Name','nlms_single_binary');
   imshow(nlms_single_binary,[]);
   title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none');
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
end

%% Determine nlmsThreshold
if(isempty(ip.Results.nlmsThreshold))
    %% Attenuate meanResponse by neighbor occupancy
    nhood_filter = [1 1 1; 1 0 1; 1 1 1];
    nhood_occupancy = imfilter(double(nlms_single_binary),nhood_filter)/8;
    % double the occupancy for accelerated attenuation
    nhood_occupancy = nhood_occupancy*2; 
    attenuatedMeanResponse = (1-nhood_occupancy).*meanResponse;
    attenuatedMeanResponse = max(attenuatedMeanResponse,0);

    nlmsThreshold = attenuatedMeanResponse;
else
    % User defined nlmsThreshold
    nlmsThreshold = ip.Results.nlmsThreshold;
end

if(ip.Results.diagnosticMode)
    other.diagFig(end+1) = figure('Name','nlmsThreshold');
    imshow(nlmsThreshold,[]);
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none');
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
end


%% Calculate high resolution maxima

% Adapt length
if(ip.Results.adaptLengthInRegime)
    % Find orientation maxima with nlmsMask only
    [maxima_highest_temp,minima_highest_temp] = interpft_extrema(a_hat,1,true,[],false);
    % Count
    n_maxima_highest_temp = ...
        size(maxima_highest_temp,1) - sum(isnan(maxima_highest_temp));
    K_high = F.K;
    K_low = max(n_maxima_highest_temp - 1,ip.Results.responseOrder);
    [K_high,K_low] = ...
        orientationSpace.diffusion.findRegimeBifurcation( ...
            a_hat,F.K, ...
            K_high,K_low, ...
            maxima_highest_temp,minima_highest_temp, ...
            [],0.1,true);
    best_derivs = orientationSpace.diffusion.orientationMaximaFirstDerivative(a_hat,F.K,maxima_highest_temp);
    best_abs_derivs = abs(best_derivs);
    best_K = repmat(F.K,size(best_derivs));
    best_maxima = maxima_highest_temp;
    maxima_working = maxima_highest_temp;
    for K=F.K:-ip.Results.K_sampling_delta:1
        s = K > K_high;
        lower_a_hat = orientationSpace.getResponseAtOrderVecHat(a_hat(:,s),F.K,K);
        [new_derivs(:,s),~,maxima_working(:,s)] = orientationSpace.diffusion.orientationMaximaFirstDerivative(lower_a_hat,K,maxima_working(:,s),[],true);
        new_abs_derivs(:,s) = abs(new_derivs(:,s));
        better(:,s) = new_abs_derivs(:,s) < best_abs_derivs(:,s);
        
        % Update better
        best_abs_derivs(better) = new_abs_derivs(better);
        best_derivs(better) = new_derivs(better);
        best_K(better) = K;
        best_maxima(better) = maxima_working(better);
    end
    
    maxima_highest_temp = best_maxima / 2;
else
    % Find orientation maxima with nlmsMask only
    maxima_highest_temp = interpft_extrema(a_hat,1,true,[],false)/2;
    best_K = repmat(F.K,size(maxima_highest_temp));
end

maxima_highest = nanTemplate(:,:,ones(size(maxima_highest_temp,1),1));
maxima_highest = shiftdim(maxima_highest,2);
for i=1:size(maxima_highest_temp,1)
    maxima_highest(i,nlmsMask) = maxima_highest_temp(i,:);
end
maxima_highest = shiftdim(maxima_highest,1);
if(ip.Results.diagnosticMode)
    other.diagFig(end+1) = figure('Name','maxima_highest(:,:,1)');
    imshow(maxima_highest(:,:,1),[]);
    caxis([0 pi]); colormap(hsv);
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none')
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
    
    other.diagFig(end+1) = figure('Name','maxima_highest');
    cm = colorcet('CBC2','N',256);
    maxima_highest_response = R_res.interpft1(maxima_highest);
    imshow(max(maxima_highest_response,[],3),[]);
    orientationSpace.rainbowOrientationQuivers(maxima_highest,maxima_highest_response,cm);
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none');
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
end

%% Perform NLMS

% Perform nlms within nlmsMask and apply nlmsThreshold for highest
% resolution maxima at lower response level
nlms_highest = nonLocalMaximaSuppressionPrecise(real(R_res.a),maxima_highest,[],[],nlmsMask);
nlms_highest_mip = max(nlms_highest,[],3);
nlms_highest_mip_binary = nlms_highest_mip > nlmsThreshold;
if(ip.Results.diagnosticMode)
    other.diagFig(end+1) = figure('Name','nlms_highest_mip');
    imshow(nlms_highest_mip,[]);
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none')
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
    
    other.diagFig(end+1) = figure('Name','nlms_highest_mip_binary');
    imshow(nlms_highest_mip_binary,[]);
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none');
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
end

% Perform nlms within nlmsMask and apply nlmsThreshold for highest
% resolution maxima and single maxima at high response level
nlmsR = nonLocalMaximaSuppressionPrecise(real(R.a),cat(3,maxima_highest,maximum_single_angle),[],[],nlmsMask);
nlmsR_mip = max(nlmsR,[],3);
nlmsR_mip_binary = nlmsR_mip > nlmsThreshold;
if(ip.Results.diagnosticMode)
    other.diagFig(end+1) = figure('Name','nlmsR_mip');
    imshow(nlmsR_mip,[]);
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none');
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
    
    other.diagFig(end+1) = figure('Name','nlmsR_mip_binary');
    imshow(nlmsR_mip_binary,[]);
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none');
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
end

% Calculate the maximum projected response at the lower response level
% using highest and lowest resolution maxima
combinedMaxima = cat(3,maxima_highest,maximum_single_angle);
combinedResR = R_res.interpft1(combinedMaxima);
combinedResR_mip = max(combinedResR,[],3);
if(ip.Results.diagnosticMode)
    other.diagFig(end+1) = figure('Name','combinedResR_mip');
    imshow(combinedResR_mip,[]);
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none');
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
end

%% Calculate consensus and determine pixels for potential improvement
consensus_binary = nlms_highest_mip_binary & nlms_single_binary;
if(ip.Results.diagnosticMode)
    other.diagFig(end+1) = figure('Name','consensus_binary');
    imshow(consensus_binary,[]);
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none');
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
end

%% Do minimal bridging

if(~ip.Results.diagnosticMode)
    [bridgedSkeletonOne, other.bridging] = ...
        bridging.minimalBridgingFromBinary( ...
        I, ...
        nlms_highest_mip_binary, ...
        consensus_binary, ...
        5, ...
        ip.Results.diagnosticMode);
else

    [bridgedSkeletonOne, other.bridging, other.diagFig(end+1:end+2)] = ...
        bridging.minimalBridgingFromBinary( ...
        I, ...
        nlms_highest_mip_binary, ...
        consensus_binary, ...
        5, ...
        ip.Results.diagnosticMode, ...
        diag_rp.BoundingBox);
    
end

%% Repeat bridging with nlmsR

if(~ip.Results.diagnosticMode)
    [bridgedSkeletonTwo, other.bridging(2)] = ...
            bridging.minimalBridgingFromBinary( ...
            I, ...
            nlmsR_mip_binary | bridgedSkeletonOne, ...
            bridgedSkeletonOne, ...
            10, ...
            ip.Results.diagnosticMode);
else

    [bridgedSkeletonTwo, other.bridging(2), other.diagFig(end+1:end+2)] = ...
        bridging.minimalBridgingFromBinary( ...
        I, ...
        nlmsR_mip_binary | bridgedSkeletonOne, ...
        bridgedSkeletonOne, ...
        10, ...
        ip.Results.diagnosticMode, ...
        diag_rp.BoundingBox);

    
end



response = combinedResR;
theta = combinedMaxima;
response_nan = isnan(response);
response(response_nan) = -Inf;
[response,theta,response_nan] = sortMatrices(response,theta,response_nan,3,'descend');
response(response_nan) = NaN;

nms = bridgedSkeletonTwo.*combinedResR_mip;
angularResponse = R.a;

if(nargout > 4)
    % Process and output NLMS / maxima from highest K
    nlms_highest_nan = isnan(nlms_highest);
    nlms_highest(nlms_highest_nan) = -Inf;
    [nlms_highest,maxima_highest,nlms_highest_nan,combinedResR(:,:,1:end-1)] = sortMatrices(nlms_highest,maxima_highest,nlms_highest_nan,combinedResR(:,:,1:end-1),3,'descend');
    nlms_highest(nlms_highest_nan) = NaN;
    
    other.nlms_highest = nlms_highest;
    other.maxima_highest = maxima_highest;
    other.nlms_highest_mip = nlms_highest_mip;
    
    % Process and output NLMS / maxima from single regime
    other.maximum_single_angle = maximum_single_angle;
    other.nlms_single = nlms_single;
    other.nlms_single_binary = nlms_single_binary;
    other.attenuatedMeanResponse = attenuatedMeanResponse;
    other.meanResponse = meanResponse;
    other.nlmsMask = nlmsMask;
    other.params = ip.Results;
    other.combinedResR = combinedResR;
    other.nlmsR = nlmsR;
    other.nlmsR_mip_binary = nlmsR_mip_binary;
    
    % format K array
    other.K_highest = nanTemplate(:,:,ones(size(best_K,1),1));
    other.K_highest = shiftdim(other.K_highest,2);
    for i=1:size(best_K,1)
        other.K_highest(i,nlmsMask) = best_K(i,:);
    end
    other.K_highest = shiftdim(other.K_highest,1);
    other.K_highest(isnan(maxima_highest)) = NaN;
end

if(ip.Results.diagnosticMode)
    other.diagFig(end+1) = figure('Name','response');
    imshow(max(response,[],3),[]);
    h = orientationSpace.rainbowOrientationQuivers(theta,response);
    set(h,'LineWidth',2);
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none');
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
    
    other.diagFig(end+1) = figure('Name','nms');
    imshow(nms,[]);    
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none')
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
    
    other.diagFig(end+1) = figure('Name','nms > 0');
    imshow(nms > 0,[]);    
    title(other.diagFig(end).Children,other.diagFig(end).Name, ...
        'Interpreter','none')
    xlim([0 diag_rp.BoundingBox(3)]+diag_rp.BoundingBox(1));
    ylim([0 diag_rp.BoundingBox(4)]+diag_rp.BoundingBox(2));
    
    iptsetpref('ImshowBorder',oldPrefImshowBorder);
end


end
steerableAdaptiveLengthOrientationSpaceDetector Perform adaptive length orientation
space segmentation

 INPUT
 I - image
     Type: 2D numeric matrix, non-empty
 order - (optional), K parameter for OrientationSpaceFilter
     Type: Numeric, scalar
     Default: 8
 sigma - (optional), scale parameter setting the radial bandpass in pixels
         central frequency of the bandpass fliter will be 1/2/pi/sigma
     Type: Numeric, scalar
     Default: 2 (pixels)

 PARAMETERS
 filter - OrientationSpaceFilter to use, overrides order and sigma
 response - OrientationSpaceResponse to use, overrides order and sigma
 adaptLengthInRegime - Search current maxima regime to find smallest slopes
     Type: logical
     Default: true
 meanThresholdMethod - function to determine threshold of mean response
     Type: char, function\_handle
     Default: @thresholdOtsu
 mask - binary mask the same size as I
     Type: logical
     Default: []
 nlmsMask - binary mask the same size as I to apply to NLMS
     Type: logical
     Default: []
 nlmsThreshold - threshold to apply to NLMS
     Type: numeric, 2D
     Default: []
 useParallelPool - true if parallel pool should be used
     Type: logical
     Default: true
 maskDilationDiskRadius - radius in pixels to dilate the mask 
     Type: numeric
     Default: 3
 maskFillHoles - fill holes in the mask
     Type: logical
     Default: false
 diagnosticMode - true if diagnostic figures should be shown
     Type: logical, scalar
     Default: false
 K\_sampling\_delta - interval to sample K when using adaptLengthInRegime
     Type: numeric, scalar
     Default: 0.1
 responseOrder - order from which to retrieve response (not maxima)
     Type: numeric, scalar
     Default: 3

 OUTPUT
 response - response values (at responseOrder) of maxima
            size(image) x M
 theta    - maxima
 nms      - non-maximum suppression like image
            same size as image
 angularResponse
          - response values at regular intervals
 other    - struct containing various variables:
 .nlms\_highest - NLMS using highest regime maxima
 .maxima\_highest - Maxima at highest regime
 .nlms\_highest\_mip - Maximum Intensity Projection of .nlms\_highest
 .maximum\_single\_angle - Angle at lowest single-orientation regime
 .nlms\_single - NLMS at maximum\_single\_angle
 .nlms\_single\_binary - binarized .nlms\_single by meanResponse
 .attenuatedMeanResponse - meanResponse attenuated by neighborhood
 .meanResponse - mean orientationSpace response
 .nlmsMask - mask used to process NLMS
 .params - input parameters
 .combinedResR - currently the same as response


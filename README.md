# Adaptive Resolution Orientation Space

## Abstract
Microscopy images of cytoskeletal, nucleoskeletal, and other filamentous
structures contain complex junctions where multiple filaments with distinct
orientations overlap, yet state-of-the-art software generally uses single
orientation analysis to segment these structures.
We describe an image analysis approach to simultaneously segment both
filamentous structures and their intersections in microscopy images,
based on analytically resolving coincident multiple orientations upon
image filtering in a manner that balances orientation resolution and
spatial localization.


## Manuscript

"Adaptive-Resolution Multi-Orientation Analysis of Complex Filamentous Network Images"

Mark Kittisopikul<sup>1,2</sup>, Amir Vahabikashi<sup>2</sup>, Takeshi Shimi<sup>2,3</sup>, Robert D. Goldman<sup>2</sup>, and Khuloud Jaqaman<sup>1,4</sup>

### Affiliations

1. Department of Biophysics, UT Southwestern Medical Center, Dallas, TX 75390
2. Department of Cell and Developmental Biology, Feinberg School of Medicine, Northwestern University, Chicago, IL 60611
3. Institute of Innovative Research, Tokyo Institute of Technology, Yokohama, Japan
4. Department of Bioinformatics, UT Southwestern Medical Center, Dallas, TX 75390


## steerableAdaptiveLengthOrientationSpaceDetector
Perform adaptive length orientation space segmentation

### INPUT
    I - image
        Type: 2D numeric matrix, non-empty
    order - (optional), K parameter for OrientationSpaceFilter
        Type: Numeric, scalar
        Default: 8
    sigma - (optional), scale parameter setting the radial bandpass in pixels
            central frequency of the bandpass fliter will be 1/2/pi/sigma
        Type: Numeric, scalar
        Default: 2 (pixels)

### PARAMETERS
    filter - OrientationSpaceFilter to use, overrides order and sigma
    response - OrientationSpaceResponse to use, overrides order and sigma
    adaptLengthInRegime - Search current maxima regime to find smallest slopes
        Type: logical
        Default: true
    meanThresholdMethod - function to determine threshold of mean response
        Type: char, function_handle
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
    K_sampling_delta - interval to sample K when using adaptLengthInRegime
        Type: numeric, scalar
        Default: 0.1
    responseOrder - order from which to retrieve response (not maxima)
        Type: numeric, scalar
        Default: 3

### OUTPUT
    response - response values (at responseOrder) of maxima
               size(image) x M
    theta    - maxima
    nms      - non-maximum suppression like image
               same size as image
    angularResponse
             - response values at regular intervals
    other    - struct containing various variables:
    .nlms_highest - NLMS using highest regime maxima
    .maxima_highest - Maxima at highest regime
    .nlms_highest_mip - Maximum Intensity Projection of .nlms_highest
    .maximum_single_angle - Angle at lowest single-orientation regime
    .nlms_single - NLMS at maximum_single_angle
    .nlms_single_binary - binarized .nlms_single by meanResponse
    .attenuatedMeanResponse - meanResponse attenuated by neighborhood
    .meanResponse - mean orientationSpace response
    .nlmsMask - mask used to process NLMS
    .params - input parameters
    .combinedResR - currently the same as response

### EXAMPLES
    demo = zeros(256);
    demo(128,:) = 1;
    demo = max(imgaussfilt(demo,2),imgaussfilt(eye(256),2));
    demo = imnoise(mat2gray(demo),'gaussian',0.1,0.01);
    [res,theta,nms] = steerableAdaptiveLengthOrientationSpaceDetector(demo);
    figure; imshow(nms,[]);
    orientationSpace.rainbowOrientationQuivers(theta,res,hsv(32));
    xlim(128+[-10 10]);
    ylim(128+[-10 10]);

![Zoom in Demonstration of Adaptive Resolution Orientation Space and NLMS Analysis](demo/demo.png)

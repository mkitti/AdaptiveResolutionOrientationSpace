    demo = zeros(256);
    demo(128,:) = 1;
    demo = max(imgaussfilt(demo,2),imgaussfilt(eye(256),2));
    demo = imnoise(mat2gray(demo),'gaussian',0.1,0.01);
    [res,theta,nms] = steerableAdaptiveLengthOrientationSpaceDetector(demo);
    figure; imshow(nms,[]);
    orientationSpace.rainbowOrientationQuivers(theta,res,hsv(32));
    xlim(128+[-10 10]);
    ylim(128+[-10 10]);
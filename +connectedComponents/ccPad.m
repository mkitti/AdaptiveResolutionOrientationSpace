function padded = ccPad(cc,padSize)
    padded = cc;
    padded.ImageSize = cc.ImageSize + padSize*2;
    ccImageSize = cc.ImageSize;
    paddedImageSize = padded.ImageSize;
    padded.PixelIdxList = cellfun(@padPixels,cc.PixelIdxList,'UniformOutput',false);
%     for i=1:cc.NumObjects
%         [r,c] = ind2sub(cc.ImageSize , cc.PixelIdxList{i});
%         r = r + padSize(1);
%         c = c + padSize(2);
%         padded.PixelIdxList{i} = sub2ind(padded.ImageSize , r, c);
%     end
    padded.Pad = padSize;
    function paddedPixels = padPixels(pixelIdx)
        [r,c] = ind2sub(ccImageSize, pixelIdx);
        r = r + padSize(1);
        c = c + padSize(2);
        paddedPixels = sub2ind(paddedImageSize , r, c);
    end
end

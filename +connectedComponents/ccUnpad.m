function padded = ccUnpad(cc,padSize)
    if(nargin < 2)
        padSize = cc.Pad;
    end
    padded = cc;
    padded.ImageSize = cc.ImageSize - padSize*2;
    ccImageSize = cc.ImageSize;
    paddedImageSize = padded.ImageSize;
    padded.PixelIdxList= cellfun(@unpadPixels,cc.PixelIdxList,'UniformOutput',false);
%     for i=1:cc.NumObjects
%         [r,c] = ind2sub(cc.ImageSize , cc.PixelIdxList{i});
%         r = r - padSize(1);
%         c = c - padSize(2);
%         f = r > 0 & ...
%             c > 0 & ...
%             r <= padded.ImageSize(1) & ...
%             c <= padded.ImageSize(2);
%         padded.PixelIdxList{i} = sub2ind(padded.ImageSize , r(f), c(f));
%     end
    padded = rmfield(padded,'Pad');
    function unPaddedPixels = unpadPixels(pixelIdx)
        [r,c] = ind2sub(ccImageSize , pixelIdx);
        r = r - padSize(1);
        c = c - padSize(2);
        f = r > 0 & ...
            c > 0 & ...
            r <= paddedImageSize(1) & ...
            c <= paddedImageSize(2);
        unPaddedPixels = sub2ind(paddedImageSize , r(f), c(f));
    end
end

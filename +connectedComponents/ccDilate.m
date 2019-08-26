function dilatedCC = ccDilate(cc,se)
    import connectedComponents.*;
    if(isa(se,'strel'))
        se = se.getnhood;
    end
    
%     cc.PixelIdxList = distributed(cc.PixelIdxList);
    
    dilatedCC = ccPad(cc,size(se));
    offsets = mat2offset(dilatedCC.ImageSize,se);
%     test = dilatedCC.PixelIdxList;
%     test = distributed(test);

    dilatedCC.PixelIdxList = cellfun(@applyOffsets,dilatedCC.PixelIdxList,'UniformOutput',false);
%     for i=1:cc.NumObjects
%         dilatedCC.PixelIdxList{i} = bsxfun(@plus,dilatedCC.PixelIdxList{i},offsets);
%         dilatedCC.PixelIdxList{i} = unique(dilatedCC.PixelIdxList{i});
%         dilatedCC.PixelIdxList{i} = dilatedCC.PixelIdxList{i}(:);
%     end
    dilatedCC = ccUnpad(dilatedCC);
%     dilatedCC.PixelIdxList = gather(dilatedCC.PixelIdxList);
    function offetPixelIdx = applyOffsets(pixelIdx)
        offetPixelIdx = bsxfun(@plus,pixelIdx,offsets);
        offetPixelIdx = unique(offetPixelIdx);
        offetPixelIdx = offetPixelIdx(:);
    end
end
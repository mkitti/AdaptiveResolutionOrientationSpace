function offsets  = mat2offset(siz,m)
    % obtain center row and column coordinates
    center = ceil(size(m)/2);
    centerInd = sub2ind(siz,center(1),center(2));
    [r,c] = find(m);
    ind = sub2ind(siz,r,c);
    offsets = ind - centerInd;
    offsets = offsets';
end

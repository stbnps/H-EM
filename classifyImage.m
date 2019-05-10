function [ res ] = classifyImage( inputImage, classifier)
%CLASSIFYIMAGE Sliding window classifier
%   The input image is divided into 9x9 patches, and every patch is
%   classified as background or cell




patchSize = 9;
[rows, columns] = size(inputImage);
stride = 1;

[xIndices, yIndices] = ndgrid( 1:stride:(columns-(patchSize - 1)) , 1:stride:(rows - (patchSize - 1)));


indices = bsxfun(@plus, (1 : patchSize)', ((1 : patchSize) - 1) * rows);

offset = (yIndices(:)-1) + (xIndices(:)-1) * rows;

% Increase this to process more patches simultaneously
pageSize = 140000;

pageCount = ceil(numel(offset) / pageSize);

pagingIndices = [0 : pageCount - 1] .* pageSize + 1;

resCell = cell(pageCount, 1);

parfor i = 1 : numel(pagingIndices)
    
    pageStart = pagingIndices(i);
    pageEnd = min([numel(offset) pageStart + pageSize - 1]);

    currentOffset = offset(pageStart : pageEnd);
    pageIndices = bsxfun(@plus, indices, reshape(currentOffset, [1 1 numel(currentOffset)]));
    
    
    % This speeds up the function by slicing multiple patches simultaneously,
    % so that the for loop doesnt iterate all pixels

    iDescriptors = inputImage(pageIndices);
    iDescriptors = reshape(iDescriptors, patchSize * patchSize, size(iDescriptors, 3))';

    minIDescriptors = min(iDescriptors, [], 2);
    maxIDescriptors = max(iDescriptors, [], 2);
    rangeIDescriptors = maxIDescriptors - minIDescriptors;
    iDescriptors = bsxfun(@minus, iDescriptors, minIDescriptors);
    descriptors = bsxfun(@rdivide, iDescriptors, rangeIDescriptors);

    [labels, probs] = predict(classifier, descriptors);
    resCell{i} = probs(:, 1);

end

res = cell2mat(resCell);


res = reshape(res, size(inputImage, 2) - (patchSize - 1), size(inputImage, 1) - (patchSize - 1))';



res = padarray(res, [(patchSize - 1)/2 (patchSize - 1)/2]);


end


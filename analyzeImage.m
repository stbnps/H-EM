function [ cellFluorescences, cellSizes, cellCoordinates ] = analyzeImage( inputImage )
%ANALYZEIMAGE Returns the fluorescence and size of the cells in a image
%   First detects every bright spot using a sliding window classifier and
%   then applies the H-EM algorithm to each spot



% First, detect bright spots

classifier = load('classifier'); % pretrained classifier
classifier = classifier.classifier;

cellMask = classifyImage(inputImage, classifier);

detectorResponseThreshold = (1 - 1e-2); 

cellMask = cellMask > detectorResponseThreshold;

cc = bwconncomp(cellMask);

cellCoordinates = regionprops(cc,'Centroid'); 

cellCoordinates = cat(1, cellCoordinates.Centroid);

cellCoordinates = fliplr(cellCoordinates);


% The cells close to the borders are discarded, since they are out of focus

borderMargin = 64;
cellCoordinates(cellCoordinates(:, 1) < borderMargin, :) = [];
cellCoordinates(cellCoordinates(:, 1) > size(inputImage, 1) - borderMargin, :) = [];
cellCoordinates(cellCoordinates(:, 2) < borderMargin, :) = [];
cellCoordinates(cellCoordinates(:, 2) > size(inputImage, 2) - borderMargin, :) = [];


cellSizes = zeros(size(cellCoordinates, 1), 1);
cellFluorescences = zeros(size(cellCoordinates, 1), 1);
cellCircularities = zeros(size(cellCoordinates, 1), 1);


% H-EM is computed for every cell locally, on 15x15 patches

patchRadius = 7;

parfor i = 1 : size(cellCoordinates, 1)

    y = round(cellCoordinates(i, 1));
    x = round(cellCoordinates(i, 2));
    patch = inputImage(y - patchRadius : y + patchRadius, x - patchRadius : x + patchRadius);

    [cellFluorescence, cellSize, cellCircularity] = computeBiomarkers(patch);
    
    

    if cellFluorescence > 0
        cellSizes(i) = cellSize;
        cellFluorescences(i) = cellFluorescence;
        cellCircularities(i) = cellCircularity;

    end
end


% Finally, remove outliers

% Large cells
cellSizeThreshold = 300;

% Bright cells
fluorescenceHighTrhreshold = 10 ^ 6;

% Dim cells
fluorescenceLowTrhreshold = 10;

% Stretched cells
circularityThreshold = 2;


cellFluorescences(cellCircularities > circularityThreshold) = [];
cellCoordinates(cellCircularities > circularityThreshold, :) = [];
cellSizes(cellCircularities > circularityThreshold) = [];
cellCircularities(cellCircularities > circularityThreshold) = [];


cellFluorescences(cellSizes > cellSizeThreshold) = [];
cellCoordinates(cellSizes > cellSizeThreshold, :) = [];
cellCircularities(cellSizes > cellSizeThreshold) = [];
cellSizes(cellSizes > cellSizeThreshold) = [];


cellCoordinates(cellFluorescences  > fluorescenceHighTrhreshold, :) = [];
cellSizes(cellFluorescences  > fluorescenceHighTrhreshold) = [];
cellCircularities(cellFluorescences  > fluorescenceHighTrhreshold) = [];
cellFluorescences(cellFluorescences  > fluorescenceHighTrhreshold) = [];



cellCoordinates(cellFluorescences  < fluorescenceLowTrhreshold, :) = [];
cellSizes(cellFluorescences  < fluorescenceLowTrhreshold) = [];
cellFluorescences(cellFluorescences  < fluorescenceLowTrhreshold) = [];


end




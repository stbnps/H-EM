function [ cellFluorescence, cellSize, circularity ] = computeBiomarkers( imPatch )
%COMPUTEBIOMARKERS Computes biomarkers for a cell in a given patch
%   Fluorescence intensity, size, and circularity are computed using the
%   H-EM algorithm

nGaussians = 1;

mu = size(imPatch) / 2;

% H-EM returns a Gausian + Uniform Mixture (GUM)
gum = hem(imPatch, nGaussians, mu, true);


[cellFluorescences, cellProbabilities, marginalProbability] = gum2fluorescences(gum, imPatch);

if size(cellFluorescences) == 1
    sigma = gum.sigma(:, :, 1);
    [~, eigenValues] = eig(sigma);
    eigenValues = diag(eigenValues);
    a = eigenValues(1);
    b = eigenValues(2);
    cellFluorescence = cellFluorescences(1);
    cellSize = 6 * sqrt(a) * 6 * sqrt(b);
    circularity = max(a,b) / min(a,b);
else
    cellFluorescence = -1;
    cellSize = -1;    
    circularity = -1;
end


end








function [ cellFluorescences, cellProbabilities, marginalProbability] = gum2fluorescences( gum, imPatch)
%GUM2FLUORESCENCES Computes the fluorescence intensity given a GUM
%   Integrates de probability of the graussian distribution

muRecovered = gum.mu;
sigmaRecovered = gum.sigma;

minY = gum.minY;
maxY = gum.maxY;
minX = gum.minX;
maxX = gum.maxX;

cellFluorescences = zeros(size(muRecovered, 1), 1);

marginalProbability = getGUMMarginal(gum);

cellProbabilities = zeros(size(imPatch));


range = 3;
for i = 1 : size(muRecovered, 1)

    s = max(diag(sigmaRecovered(:, :, i)));

    lX = ceil(muRecovered(i, 2) - range * s);
    if lX < 1
        lX = 1;
    end

    rX = floor(muRecovered(i, 2) + range * s);
    if rX > maxX
        rX = maxX;
    end

    uY = ceil(muRecovered(i, 1) - range * s);
    if uY < 1
        uY = 1;
    end

    dY = floor(muRecovered(i, 1) + range * s);
    if dY > maxY
        dY = maxY;
    end

    [x, y] = meshgrid(lX : rX,  uY : dY);

    data_i = [y(:) x(:)];
    closeIndices = sub2ind([maxY, maxX], y(:), x(:));

    try
        f = mvnpdf(data_i, muRecovered(i, :), sigmaRecovered(:, :, i));
    catch
        continue;
    end


    cellProbability = gum.pi(i) .* f ./ marginalProbability(closeIndices);

    cellProbabilities(closeIndices) = cellProbabilities(closeIndices) + cellProbability(:);

    cell = cellProbability .* imPatch(closeIndices);

    cellFluorescences(i) = sum(cell(:));

end


end


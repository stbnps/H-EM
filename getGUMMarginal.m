function [ marginal ] = getGUMMarginal( gum )
%GETGUMMARGINAL Computes the marginal probability of a GUM
%   -

minY = gum.minY;
maxY = gum.maxY;
minX = gum.minX;
maxX = gum.maxX;

uniformP = 1 / ((maxY - minY + 1) * (maxX - minX + 1));
uniformP = repmat(uniformP, (maxY - minY + 1) * (maxX - minX + 1), 1);
mu = gum.mu;
sigma = gum.sigma;

marginal = zeros(maxY - minY + 1, maxX - minX + 1);

range = 3;

for i = 1 : gum.gaussianDistributionCount
    
    s = max(diag(sigma(:, :, i)));

    lX = ceil(mu(i, 2) - range * s);
    if lX < 1
        lX = 1;
    end

    rX = floor(mu(i, 2) + range * s);
    if rX > maxX
        rX = maxX;
    end

    uY = ceil(mu(i, 1) - range * s);
    if uY < 1
        uY = 1;
    end

    dY = floor(mu(i, 1) + range * s);
    if dY > maxY
        dY = maxY;
    end

    [x, y] = meshgrid(lX : rX,  uY : dY);
    closeIndices = sub2ind([maxY, maxX], y(:), x(:));
    
    f = zeros(maxY - minY + 1, maxX - minX + 1);

    try
        f_i = mvnpdf([y(:), x(:)], mu(i, :), sigma(:, :, i));
    catch
        continue;
    end

    f_i = f_i .* gum.pi(i);
    f(closeIndices) = f_i(:);
    
    marginal(closeIndices) = marginal(closeIndices) + f(closeIndices);

end

f = gum.pi(gum.gaussianDistributionCount + 1) .* uniformP;
f = reshape(f, maxY - minY + 1, maxX - minX + 1);
marginal = marginal + f;

end


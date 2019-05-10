function [gum] = hem( inputImage, gaussianDistributionCount, mu, trimming)
%HEM Fits H-EM to an image patch
%   Performs expectation maximization over a mixture model using an uniform
%   and gaussian distributions

[x, y] = meshgrid(1:size(inputImage,2), 1:size(inputImage, 1));

% In this function "event" means a increase of 1 in the intensity value of
% a given pixel
eventLocations = [y(:) x(:)];
eventCounts = inputImage(1:end);



tolf = 1e-6;
maxIterations = 32;

minY = min(eventLocations(:, 1));
maxY = max(eventLocations(:, 1));
minX = min(eventLocations(:, 2));
maxX = max(eventLocations(:, 2));
    
dataCount = size(eventLocations, 1);    
eventCountsSum = sum(eventCounts(:));

distributionCount = gaussianDistributionCount + 1;
pi = ones(distributionCount, 1) ./ distributionCount; 


if nargin < 3
    mu = rand(gaussianDistributionCount, 2) .* repmat([maxY - minY, maxX - minX], gaussianDistributionCount, 1) + repmat([minY, minX], gaussianDistributionCount, 1);
end

sigma = repmat(2 * eye(2), 1, 1, gaussianDistributionCount);

finished = false;

likelihood = -inf;

iterationCount = 0;

expectation = cell(distributionCount, 1);

processedIndices = cell(gaussianDistributionCount, 1);


% Perform expectation maximization until convergence
while ~finished

    previousLikelihood = likelihood;

    deletedElementCount = 0;
    if trimming
        deletedIndices = zeros(gaussianDistributionCount, 1);
    end
    
    % The gaussian distributions are only evaluated on +-4 sigma
    range = 4;
    marginalSum = zeros(1, size(eventLocations, 1));
    for i = 1 : gaussianDistributionCount
        

        
        s = sqrt(max(diag(sigma(:, :, i))));
        
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

        data_i = [y(:) x(:)];
        closeIndices = sub2ind([maxY, maxX], y(:), x(:));
        
        if trimming
        
            if numel(closeIndices(:)) <= 1
                deletedElementCount = deletedElementCount + 1;
                deletedIndices(deletedElementCount) = i;
                continue;
            end

            if sigma(1, 1, i) < 1.01/3 || sigma(2, 2, i) < 1.01/3
                deletedElementCount = deletedElementCount + 1;
                deletedIndices(deletedElementCount) = i;
                continue;
            end

            sigmaLimit = 24;
            if sigma(1, 1, i) > sigmaLimit || sigma(2, 2, i) > sigmaLimit
                deletedElementCount = deletedElementCount + 1;
                deletedIndices(deletedElementCount) = i;
                continue;
            end  
        
        end

        processedIndices{i} = closeIndices;

        
        
        try
            expectation{i} = pi(i) .* mvnpdf(data_i, mu(i, :), sigma(:, :, i))';
        catch
            deletedElementCount = deletedElementCount + 1;
            deletedIndices(deletedElementCount) = i;
            continue;
        end
        
        marginalSum(closeIndices) = marginalSum(closeIndices) + expectation{i};


    end
    
    % delete unwanted data
    if trimming && deletedElementCount > 0
        deletedIndices = deletedIndices(1 : deletedElementCount);
        pi(deletedIndices) = [];
        mu(deletedIndices, :) = [];
        sigma(:, :, deletedIndices) = [];
        expectation(deletedIndices) = [];
        processedIndices(deletedIndices) = []; 
        distributionCount = distributionCount - deletedElementCount;
        gaussianDistributionCount = gaussianDistributionCount - deletedElementCount;
    end
    
    % Calculate the uniform part
    uniformComponent = 1 / ((1 + maxX - minX)  * (1 + maxY - minY));
    expectation{distributionCount} = pi(distributionCount) .* uniformComponent;

    marginalSum = marginalSum + pi(distributionCount) .* uniformComponent;



    likelihood = sum(log(marginalSum) .* eventCounts, 2);


    for i = 1 : size(expectation, 1) - 1
        
        expectation{i} = expectation{i}./ marginalSum(processedIndices{i});

        
    end
    
	expectation{distributionCount} = expectation{distributionCount} ./ marginalSum(1, :);


    % Maximization:
    for i = 1 : gaussianDistributionCount
        

        
        % Cache variables by hand
        
        expectationTimesEventcounts = expectation{i} .* eventCounts(processedIndices{i});

        expectationTimesEventcountsTimesData = bsxfun(@times, eventLocations(processedIndices{i}, :), expectationTimesEventcounts');
        expectationTimesEventcountsTimesDataSum = sum(expectationTimesEventcountsTimesData);
        expectationTimesEventcountsSum = sum(expectationTimesEventcounts);
        dataMinusMu = bsxfun(@minus, eventLocations(processedIndices{i}, :), mu(i, :));

        pi(i) = expectationTimesEventcountsSum / eventCountsSum;

        
        mu(i, :) = expectationTimesEventcountsTimesDataSum ./ expectationTimesEventcountsSum;
        
        sigma_g = bsxfun(@times, dataMinusMu, expectationTimesEventcounts')' * dataMinusMu ./ expectationTimesEventcountsSum;


        sigma(:, :, i) = sigma_g;




    end
    
    pi(distributionCount) = sum(expectation{distributionCount} .* eventCounts) / eventCountsSum;



    iterationCount = iterationCount + 1;
    
    
    if likelihood - previousLikelihood < tolf * abs(previousLikelihood) || iterationCount > maxIterations 

        finished = true;

    end




end    


gum.gaussianDistributionCount = gaussianDistributionCount;
gum.pi = pi;
gum.mu = mu;
gum.sigma = sigma;
gum.minY = minY;
gum.maxY = maxY;
gum.minX = minX;
gum.maxX = maxX;

end


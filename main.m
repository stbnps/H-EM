
dataPaths = {'test_data/Granulocyte'; 'test_data/Monocyte'; 'test_data/Lymphocyte'};

fluorescences = [];
sizes = [];
labels = [];

% Perform analysis
for i = 1 : size(dataPaths, 1)

    dataPath = dataPaths{i};
    [fls, szs] = analyzeFolder(dataPath);
    
    fluorescences = [fluorescences; fls];
    sizes = [sizes; szs];
    labels = [labels; i * ones(size(fls, 1), 1)];   

end


% Plot data
for i = 1 : 3
    
    fl = fluorescences(i == labels);
    sz = sizes(i == labels);
    
    % Perform density estimates so that data points on populated areas
    % are displayed lager
    scatter_sz = ksdensity([sz, log10(fl)], [sz, log10(fl)]) * 16;
    
    % Set a different color to each cell type
    scatter_cl = [0, 0, 0];
    scatter_cl(i) = 1;    
    scatter_cl = repmat(scatter_cl, [size(fl, 1), 1]);  
    
    scatter(sz, log10(fl), scatter_sz, scatter_cl, 'filled');
    hold on;

    
end

title('Estimated fluorescence intensity and size');
xlabel('Size');
ylabel('CD45 - Alexa Fluor 532');
legend('Granulocytes','Monocytes','Lymphocytes');



function [fluorescences, sizes] = analyzeFolder(dataPath)
%ANALYZEFOLDER Analyzes all fields of view in a folder
%   Each folder may have multiple fields of view, and each field of view
%   multiple images


fileType = 'tif';

imageFolders = dir([dataPath '/FOV*']);  

nFolders = length(imageFolders);

images = cell(nFolders, 1);
folders = cell(nFolders, 1);

fluorescences = [];
sizes = [];

for i = 1 : nFolders

    curentFolder = imageFolders(i).name;
    
    folders{i} = fullfile(dataPath, curentFolder);
    
    files = dir([fullfile(dataPath, curentFolder) '/*.' fileType]);  

    nFiles = length(files);
    for j = 1 : nFiles
        currentFilename = files(j).name;
        currentImage = imread(fullfile(dataPath, curentFolder, currentFilename));
        
        % There may be multiple images of the same FOV
        % Add them to increase SNR
        if numel(images{i} > 0)
            images{i} = images{i} + double(currentImage);
        else 
            images{i} = double(currentImage);
        end
    end

end


for i = 1 : nFolders
    inputImage = images{i};

    [fl, sz, loc] = analyzeImage( inputImage );

    fluorescences = [fluorescences; fl];
    sizes = [sizes; sz];
end


end


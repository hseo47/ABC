%% 1. Specify Input Image and Prepare Output Filename
filename = 'filename.jpg'; % Replace with your filename
[~, baseFilename, ~] = fileparts(filename);
outputTxtFile = [baseFilename, '_coverage_data.txt'];

%% 2. Read/Preprocess
I = imread(filename);
if size(I, 3) == 3
    I_gray = rgb2gray(I);
else
    I_gray = I;
end

% Cropping the image to remove information bar
cropHeight = 0.90; 
I_gray = I_gray(1:round(size(I_gray, 1) * cropHeight), :);

% Enhance contrast (adaptive histogram equalization)
I_contrast = adapthisteq(I_gray, 'ClipLimit', 0.02, 'Distribution', 'rayleigh');

% Noise reduction
I_smooth = imgaussfilt(I_contrast, 2);

% Background subtraction
bg = imopen(I_smooth, strel('disk', 30));
I_sub = imsubtract(I_smooth, bg);

% Remove faint noise
thresholdValue = graythresh(I_sub) * 0.8; 
BW = imbinarize(I_sub, thresholdValue);

% Morphological cleanup
BW = bwareaopen(BW, 100);
BW = imfill(BW, 'holes');

%% 3) Refine Aggregated Bacteria Detection Using Watershed
D = -bwdist(~BW);    % Distance transform
D(~BW) = -Inf;       % Set background to -Inf
L = watershed(D);    % Apply watershed transform
BW(L == 0) = 0;      % Remove watershed boundaries

% Filter connected components by size/shape
CC = bwconncomp(BW);
stats = regionprops(CC, 'Area', 'Eccentricity', 'Solidity');

% Define thresholds
minArea = 150;          % Minimum area of bacteria
maxArea = 3000;         % Maximum area (for large aggregates)
maxEccentricity = 0.9;  
minSolidity      = 0.75;

% Keep objects that meet all thresholds
validIdx = find([stats.Area] >= minArea & ...
                [stats.Area] <= maxArea & ...
                [stats.Eccentricity] <= maxEccentricity & ...
                [stats.Solidity] >= minSolidity);
BW_refined = ismember(labelmatrix(CC), validIdx);

%% 4) Final Counting
CC_final = bwconncomp(BW_refined);
numBacteria = CC_final.NumObjects;

%% 5) Coverage Calculations
totalPixels     = numel(BW_refined);  % total pixels in cropped image
bacteriaPixels  = sum(BW_refined(:)); % total "bacteria" pixels
coveragePercent = (bacteriaPixels / totalPixels) * 100;

% Classify each connected component as single vs. aggregate
singleAggThreshold = 800;  
statsFinal = regionprops(CC_final, 'Area', 'PixelIdxList');

singlePixelCount = 0;
aggPixelCount    = 0;

for i = 1:length(statsFinal)
    if statsFinal(i).Area < singleAggThreshold
        singlePixelCount = singlePixelCount + statsFinal(i).Area;
    else
        aggPixelCount = aggPixelCount + statsFinal(i).Area;
    end
end

singleCoveragePercent = (singlePixelCount / totalPixels) * 100;
aggCoveragePercent    = (aggPixelCount    / totalPixels) * 100;

%% 6) Print Results
fprintf('Image Analyzed: %s\n', filename);
fprintf('Number of detected bacteria: %d\n', numBacteria);
fprintf('Total surface coverage by bacteria: %.2f%%\n', coveragePercent);
fprintf('Single-bacteria coverage:     %.2f%% of total image\n', singleCoveragePercent);
fprintf('Aggregated-bacteria coverage: %.2f%% of total image\n', aggCoveragePercent);

if bacteriaPixels > 0
    fracSingle = singlePixelCount / bacteriaPixels * 100;
    fracAgg    = aggPixelCount    / bacteriaPixels * 100;
    fprintf('(Within bacterial area: %.2f%% single, %.2f%% aggregate)\n', fracSingle, fracAgg);
end

%% 7) Save Results to a Text File
fileID = fopen(outputTxtFile, 'w');
fprintf(fileID, 'Image analyzed: %s\n', filename);
fprintf(fileID, 'Number of detected bacteria: %d\n', numBacteria);
fprintf(fileID, 'Total surface coverage by bacteria: %.2f%%\n', coveragePercent);
fprintf(fileID, 'Single-bacteria coverage:     %.2f%% of total image\n', singleCoveragePercent);
fprintf(fileID, 'Aggregated-bacteria coverage: %.2f%% of total image\n', aggCoveragePercent);
if bacteriaPixels > 0
    fprintf(fileID, '(Within bacterial area: %.2f%% single, %.2f%% aggregate)\n', fracSingle, fracAgg);
end
fclose(fileID);

fprintf('\nCoverage data saved to: %s\n', outputTxtFile);

%% 8) Visualization
figure;
set(gcf, 'Position', [100, 100, 800, 1200]);
set(0, 'DefaultAxesFontName', 'Helvetica', 'DefaultAxesFontSize', 14);

subplot(3, 2, 1);
imshow(I_gray);
title('1) Cropped Grayscale');

subplot(3, 2, 2);
imshow(I_contrast, []);
title('2) Enhanced Contrast');

subplot(3, 2, 3);
imshow(I_sub, []);
title('3) Background Subtracted');

subplot(3, 2, 4);
imshow(BW);
title('4) Binary Mask w/ Watershed');

subplot(3, 2, 5);
imshow(BW_refined);
title('5) Refined Mask');

baseColors = [hex2dec('A2')/255, hex2dec('D2')/255, hex2dec('DF')/255; ...
              hex2dec('F6')/255, hex2dec('EF')/255, hex2dec('BD')/255; ...
              hex2dec('E4')/255, hex2dec('C0')/255, hex2dec('87')/255; ...
              hex2dec('BC')/255, hex2dec('7C')/255, hex2dec('7C')/255];

L = labelmatrix(CC_final);
numLabels = max(L(:));

numBaseColors = size(baseColors, 1);
if numLabels > numBaseColors
    customMap = repmat(baseColors, ceil(numLabels/numBaseColors), 1);
    customMap = customMap(1:numLabels, :);
else
    customMap = baseColors;
end

coloredLabels = label2rgb(L, customMap, 'k', 'noshuffle');

subplot(3, 2, 6);
imshow(coloredLabels);
title(['6) Final Labeled Bacteria, Count = ', num2str(numBacteria)]);

set(gcf, 'DefaultAxesTitleFontWeight', 'normal');
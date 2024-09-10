clear
clc
close all

% - Parse output .csv files each containing measurement of single cell 
% to get meaurement at each progression bin (0, 0.1, 0.2, ... 1) for each 
% cell as a single trace
%
%
%
input_dir = '/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/2022-2023/20230111 53BP1 concentration/progression analysis csv/';

addpath(input_dir)
cd(input_dir)

output_dir = '/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/2024/20240827_eYFP_constructs_reanalyzed/output/';

% Get list of all CSV files in the input folder
filePattern = fullfile(input_dir, '*.csv');
csvFiles = dir(filePattern);

% Determine the number of files
n = length(csvFiles);

% Define progression bins (0 to 10, rounded to the nearest integer)
progressionBins = 0:1:10;
m = length(progressionBins);

% Initialize the output matrix (m by n)
result = NaN(m, n); % Use NaN to handle missing data if any

% Loop through each file
for i = 1:n
    % Load the data from the CSV file
    filename = fullfile(input_dir, csvFiles(i).name);
    data = readtable(filename);
    
    % Extract Progression and Intensity columns
    progression = data.Progression;
    intensity = data.AverageIntensity;
    
    % Round the Progression to the nearest 10th
    roundedProgression = ceil(progression./10);
    
    % Calculate the mean Intensity for each bin
    for j = 1:m
        % Find indices where the rounded progression equals the current bin
        binIndices = (roundedProgression == progressionBins(j));
        
        % Calculate mean Intensity for this bin
        if any(binIndices)
            result(j, i) = mean(intensity(binIndices));
        end
    end
end

% Display the result
disp('Output Matrix (Mean Intensity for each 10th Progression Bin):');
disp(result);

% Optionally, you can save the output matrix to a CSV file
outputFile = fullfile(output_dir, 'BP1_output_matrix.csv');
writematrix(result, outputFile);

disp(['Output matrix saved to ', outputFile]);

% calculate mean and std from data_matrix
avg_std = table(mean(result,2,'omitnan'), std(result,0,2,'omitnan'),'VariableNames',{'avg','std'});
writetable(avg_std,[output_dir 'BP1_avg_std.csv'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



normResult = NaN(m,n);
% normalized matrix
for ii = 1:n
    curr_column = result(:,ii);
    nonnanValues = curr_column(~isnan(curr_column));
    norm_factor = nonnanValues(1);
    normResult(:,ii) = result(:,ii)./norm_factor;
end

% Optionally, you can save the output matrix to a CSV file
normOutputFile = fullfile(output_dir, 'BP1_norm_output_matrix.csv');
writematrix(normResult, outputFile);

disp(['Output matrix saved to ', normOutputFile]);

% calculate mean and std from data_matrix
norm_avg_std = table(mean(normResult,2,'omitnan'), std(normResult,0,2,'omitnan'),'VariableNames',{'avg','std'});
writetable(norm_avg_std,[output_dir 'BP1_norm_avg_std.csv'])



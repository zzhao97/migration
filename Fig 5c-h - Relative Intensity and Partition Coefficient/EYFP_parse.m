% Confined Migration Paper 2024
%
%
% - Parse output .csv files from EYFP_dilute_full.m (or EYFP_dilute_back.m)
% and EYFP_full.m (or EYFP_back.m) to get meaurement at each progression
% bin (0, 0.1, 0.2, ... 1) for each cell as a single trace
%
%
clear
clc
close all

input_dir = '/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/2024/20240827_eYFP_constructs_reanalyzed/input/';
addpath(input_dir)
cd(input_dir)

output_dir = '/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/2024/20240827_eYFP_constructs_reanalyzed/output/';

% load the Excel spreadsheet
file_back = 'SART1_back_output_data.csv';
file_full = 'SART1_full_output_data.csv';
data_back = readtable(file_back);
data_full = readtable(file_full);

% calculate Progression (binned to nearest 0.1) as a new variable
prog_times_10 = 10 * (data_full.AreaFull - data_back.AreaBack) ./ data_full.AreaFull;
data_back.Progression = ceil(prog_times_10);
data_back.RelInt = (data_back.MaterialBack ./ data_back.AreaBack)./(data_full.MaterialFull ./ data_full.AreaFull);

% get unique filenames and Progression bins
uniqueFilenames = unique(data_back.Filename);
progressionBins = 0:1:10;

% initialize the result table
result = table();

% loop through each unique filename
for i = 1:length(uniqueFilenames)
    % filter data for the current filename
    currentFilename = uniqueFilenames{i};
    filteredData = data_back(strcmp(data_back.Filename, currentFilename), :);
    
    for j = 1:length(progressionBins)
        % filter data for the current progression bin
        currentBin = progressionBins(j);
        binData = filteredData(filteredData.Progression == currentBin, :);
        
        % calculate the average CoV and MaterialBack
        if ~isempty(binData)
            relInt = mean(binData.RelInt);
        else
            CoV_avg = NaN;
            relInt = NaN;
        end
        
        % append to the result table
        result = [result; {currentFilename, currentBin, relInt}];
    end
end

% assign column names to the result table
result.Properties.VariableNames = {'Filename_new', 'Progression', 'RelativeIntensity'};

% display final result table 
%result

% write the result to a new CSV file
%writetable(result, [output_dir 'SRRM1_relint_parsed.csv']);

% 
%
%
%
%
%
%
%

% organize results to matrix where each column is a replicate
m = length(progressionBins);  
n = length(uniqueFilenames); 
data_matrix = NaN(m, n);  

% populate the matrix with results
for i = 1:n
    current_filename = uniqueFilenames{i};
    
    % extract rows corresponding to the current filename
    filename_data = result(strcmp(result.Filename_new, current_filename), :);
    
    % loop through each progression bin and fill the matrix
    for j = 1:m
        current_bin = progressionBins(j);
        
        % find the row for the current progression bin
        row_idx = find(filename_data.Progression == current_bin);
        
        if ~isempty(row_idx)
            data_matrix(j, i) = filename_data.RelativeIntensity(row_idx);  % Fill the matrix with CoV_avg
        end
    end
end

% display the resulting matrix (rows = Progression bins, columns = Filenames)
disp(data_matrix);

% calculate mean and std from data_matrix
avg_std = table(mean(data_matrix,2,'omitnan'), std(data_matrix,0,2,'omitnan'),'VariableNames',{'avg','std'});
%writetable(avg_std,[output_dir 'SRRM1_avg_std.csv'])





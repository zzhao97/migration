% Confined migration paper JZ 2024
%
%
% - Input files are .csv files from H2B_half_nucleoli_excluding.m and 
% containing the coefficient of variation of 
% - Progression = (whole nucleus area - trailing area)/whole nucleus area
% - Progression is then binned to tenth percentage and average perinucleolar 
% heterochromatin intensity is averaged at each progression bin
% - Output table has values for individual cells at each progression bin in
% a column and a final table of average of all cells as well as SD
%
% 
clear; 
clc; 
close all;

input_directory='/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/2022-2023/20231128_FrontvsBack_H2B_Heterogeneity_Progression/CoV calculation/';
addpath(input_directory)
cd(input_directory)

half_files = dir([input_directory, 'H2B_CoV_nucleoli_excluding_front_*.csv']);

full_files = dir([input_directory, 'H2B_CoV_nucleoli_excluding_cell*.csv']);

addpath('/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/matlab/confined migration paper/H2B rear-front scripts/')

fnames_half = {};
numfiles_half = length(half_files);
[fnames_half{1:numfiles_half,1}] = half_files(:).name;

fnames_whole = {};
numfiles_whole = length(full_files);
[fnames_whole{1:numfiles_whole,1}] = full_files(:).name;

concatProgression = [];
concatCoV = [];
concatFilename = [];

concatCellAvgProg = [];
concatCellAvgCoV = [];
concatUniqueProgName = [];

for fnum=1:numfiles_half

    half_table = readtable(fnames_half{fnum,1},'Delimiter',',');
    whole_table = readtable(fnames_whole{fnum,1},'Delimiter',',');
    
    area_half = half_table.AreaFull;
    area_whole = whole_table.AreaFull;
    
    if ~isempty(strfind(fnames_half{fnum,1},'back'))
        pctProg = (area_whole - area_half) ./ area_whole .* 100;
    elseif ~isempty(strfind(fnames_half{fnum,1},'front'))
        pctProg = area_half ./ area_whole .* 100;
    end

    progressBinnedArray = round(pctProg,-1);
 
    cov = half_table.CoV_out;
    
    concatProgression = [concatProgression; progressBinnedArray];
    concatCoV = [concatCoV; cov];

    filenameReps = cellstr(repmat(fnames_whole{fnum},[height(whole_table),1]));
    concatFilename = [concatFilename; filenameReps];

    % individual cell-average
    [uniqueVals, ~, idx] = unique(progressBinnedArray);
    sumValues = accumarray(idx, cov);
    countValues = accumarray(idx, 1);

    % Compute the average at each progression bin
    averageValues = sumValues ./ countValues;

    % store progression and average intensity pair
    concatCellAvgProg = [concatCellAvgProg; uniqueVals];
    concatCellAvgCoV = [concatCellAvgCoV; averageValues];
    uniqueProgName = cellstr(repmat(fnames_whole{fnum},[length(uniqueVals),1]));
    concatUniqueProgName = [concatUniqueProgName; uniqueProgName];

end

allProgressionTable = table(concatFilename, concatProgression, concatCoV, 'VariableNames', {'Filename', 'Progression', 'AllCoV'})
writetable(allProgressionTable, [input_directory, 'front_individual_data.csv'])

uniqueProgressionTable = table(concatUniqueProgName, concatCellAvgProg, concatCellAvgCoV, 'VariableNames', {'Filename', 'Progression', 'AverageCoV'})
writetable(uniqueProgressionTable, [input_directory, 'front_prog_avg_individual_data.csv'])

%% find the mean and std. dev. of cov at each value in the progression_grid array
progression_grid = linspace(0,100,11);

for i=1:length(progression_grid)
    nonzeroCoV = nonzeros(concatCoV(concatProgression==progression_grid(i)));
    mean_AvgCoV(i,1) = mean(nonzeroCoV,'omitnan');
    std_AvgCoV(i,1) = std(nonzeroCoV,0,'omitnan');
end

progression_grid = progression_grid';

not_nan = mean_AvgCoV(~isnan(mean_AvgCoV));

norm_mean_CoV = mean_AvgCoV./(not_nan(1));
norm_sd_CoV = std_AvgCoV./(not_nan(1));

Avg_CoV_Data=table(progression_grid, mean_AvgCoV, std_AvgCoV, norm_mean_CoV, norm_sd_CoV, 'VariableNames',{'Progression','CoVMean','CoVSTD', 'NormAvg', 'NormSTD'})
writetable(Avg_CoV_Data, [input_directory, 'front_all_avg_data.csv'])

% Confined migration paper JZ 2024
%
%
% - Input files are .csv files from perinucleolarIntensity.m and backArea.m
% containing area (um^2) of trailing half and area of whole nucleus at each
% time point
% - Progression = (whole nucleus area - trailing area)/whole nucleus area
% - Progression is then binned to tenth percentage and average perinucleolar 
% heterochromatin intensity is averaged at each progression bin
% - Output table has values for individual cells at each progression bin in 
% a column, and another output table of average of all cells as well as SD
%
% 
clear; 
clc; 
close all;

back_directory=dir('/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/2024/20240630_NPM1_perinucleolar_chromatin/backAreaOutputs/*ctrl*.csv');
whole_directory=dir('/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/2024/20240630_NPM1_perinucleolar_chromatin/perinucleolarIntensityOutputs/*ctrl*.csv');

addpath('/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/2024/20240630_NPM1_perinucleolar_chromatin/backAreaOutputs/')
addpath('/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/2024/20240630_NPM1_perinucleolar_chromatin/perinucleolarIntensityOutputs/')

output_dir = '/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/2024/20240630_NPM1_perinucleolar_chromatin/progressionAvgOutputs/';
cd(output_dir)

fnames_back = {};
numfiles_back = length(back_directory);
[fnames_back{1:numfiles_back,1}] = back_directory(:).name;

fnames_whole = {};
numfiles_whole = length(whole_directory);
[fnames_whole{1:numfiles_whole,1}] = whole_directory(:).name;

concatProgression = [];
concatAvgSurfInt = [];
concatFilename = [];

concatCellAvgProg = [];
concatCellAvgPNH = [];
concatUniqueProgName = [];

for fnum=1:numfiles_back

    back_table = readtable(fnames_back{fnum,1},'Delimiter',',');
    whole_table = readtable(fnames_whole{fnum,1},'Delimiter',',');
    
    area_back = back_table.BackArea;
    area_whole = whole_table.NucArea;
    
    pctProg = (area_whole - area_back) ./ area_whole .* 100;
    progressBinnedArray = round(pctProg,-1);
 
    avg_surf_int = whole_table.MeanSurfInt;
    
    concatProgression = [concatProgression; progressBinnedArray];
    concatAvgSurfInt = [concatAvgSurfInt; avg_surf_int];

    filenameReps = cellstr(repmat(fnames_whole{fnum},[height(whole_table),1]));
    concatFilename = [concatFilename; filenameReps];

    % individual cell-average
    [uniqueVals, ~, idx] = unique(progressBinnedArray);
    sumValues = accumarray(idx, avg_surf_int);
    countValues = accumarray(idx, 1);

    % Compute the average at each progression bin
    averageValues = sumValues ./ countValues;

    % store progression and average intensity pair
    concatCellAvgProg = [concatCellAvgProg; uniqueVals];
    concatCellAvgPNH = [concatCellAvgPNH; averageValues];
    uniqueProgName = cellstr(repmat(fnames_whole{fnum},[length(uniqueVals),1]));
    concatUniqueProgName = [concatUniqueProgName; uniqueProgName];

end

allProgressionTable = table(concatFilename, concatProgression, concatAvgSurfInt, 'VariableNames', {'Filename', 'Progression', 'AveragePNHInt'})
writetable(allProgressionTable, [output_dir, 'ctrl_individual_data.csv'])

uniqueProgressionTable = table(concatUniqueProgName, concatCellAvgProg, concatCellAvgPNH, 'VariableNames', {'Filename', 'Progression', 'AveragePNHInt'})
writetable(uniqueProgressionTable, [output_dir, 'ctrl_prog_avg_individual_data.csv'])

%% Population-average
progression_grid = linspace(0,100,11);

for i=1:length(progression_grid)
    nonzeroPNH = nonzeros(concatAvgSurfInt(concatProgression==progression_grid(i)));
    mean_AvgSurfInt(i,1) = mean(nonzeroPNH,'all');
    std_AvgSurfInt(i,1) = std(nonzeroPNH,0,'all');
end

progression_grid = progression_grid';

norm_mean_SurfInt = mean_AvgSurfInt./(mean_AvgSurfInt(1));
norm_sem = std_AvgSurfInt./(mean_AvgSurfInt(1));

Avg_PNH_Int_Data=table(progression_grid, mean_AvgSurfInt, std_AvgSurfInt, norm_mean_SurfInt, norm_sem, 'VariableNames',{'Progression','AvgSurfIntMean','AvgSurfIntSTD', 'NormAvg', 'NormSTD'})
writetable(Avg_PNH_Int_Data, [output_dir, 'ctrl_all_avg_data.csv'])

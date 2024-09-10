clear; 
clc; 
close all;

input_directory = '/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/2024/20240228_BRD4FL-miRFP_migration_pilot/single crops/';
addpath(input_directory)
cd(input_directory)

addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/MDA MB 231 confined migration features/EYFP rear-whole distribution scripts/')

back_filename = 'BRD4_back_output_data_20240306.csv';
full_filename = 'BRD4_full_output_data_20240306.csv';

back_table = readtable(back_filename);
full_table = readtable(full_filename);

rel_int = (back_table.MaterialFull./back_table.AreaFull)./(full_table.MaterialFull./full_table.AreaFull);

progression = (full_table.AreaFull - back_table.AreaFull)./full_table.AreaFull*100;
progression_binned = round((floor(progression/100) + floor((progression/100-floor(progression/100))/0.05) * 0.05) * 100);

progression_grid = linspace(0,50,11);

% find the mean and std. dev. of enrichment at each value in the
% progression_grid array
for i=1:length(progression_grid)
    mean_rel_int(i,1) = median(rel_int(progression_binned==progression_grid(i)),'omitnan');
    sem_rel_int(i,1) = std(rel_int(progression_binned==progression_grid(i)),0,'omitnan');
end

progression_grid = progression_grid';
progression_grid = progression_grid/100;

mean_rel_int= mean_rel_int./mean_rel_int(1);

Tdata=table(progression_grid, mean_rel_int, sem_rel_int,'VariableNames',{'Progression','RelIntMean','RelIntSEM'}) %generate output table
writetable(Tdata,[input_directory 'BRD4_ProgressionBinned_output_data.csv'])


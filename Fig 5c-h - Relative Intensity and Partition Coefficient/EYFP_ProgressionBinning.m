clear; 
clc; 
close all;

directory='/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/2022-2023/20231210_EYFP_H2B-miRFP_overnight/';

addpath(directory)
addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/MDA MB 231 confined migration features/EYFP rear-whole distribution scripts/')
cd(directory)

back_filename = 'EYFP_only_back_output_data.csv';
full_filename = 'EYFP_only_full_output_data.csv';

back_table = readtable(back_filename);
full_table = readtable(full_filename);

rel_int = (back_table.MaterialBack./back_table.AreaBack)./(full_table.MaterialFull./full_table.AreaFull);
%PC = back_table.PC;

progression = (full_table.AreaFull - back_table.AreaBack)./full_table.AreaFull;
%progression_binned = floor(progression) + floor((progression-floor(progression))/0.025) * 0.025;
progression_binned = round(progression,1);

progression_grid = linspace(0,1,11);

% find the mean and std. dev. of enrichment at each value in the
% progression_grid array
for i=1:length(progression_grid)
    mean_rel_int(i,1) = median(rel_int(progression_binned==progression_grid(i)),'omitnan');
    std_rel_int(i,1) = std(rel_int(progression_binned==progression_grid(i)),0,'omitnan');
    %mean_PC(i,1) = mean(PC(progression_binned==progression_grid(i)),'omitnan');
    %SE_PC(i,1) = std(PC(progression_binned==progression_grid(i)),0,'omitnan')./sqrt(length(PC(progression_binned==progression_grid(i))));
end

progression_grid = progression_grid';
Tdata=table(progression_grid, mean_rel_int, std_rel_int,'VariableNames',{'Progression','RelIntMean','RelIntSE'}) %generate output table

%Tdata=table(progression_grid, mean_PC, SE_PC,'VariableNames',{'Progression','PC_Mean','PC_SE'}) %generate output table
%writetable(Tdata,[directory 'SRRM1_ProgressionBinned_PC.csv'])


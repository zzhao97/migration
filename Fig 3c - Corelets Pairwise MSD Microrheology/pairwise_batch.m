% 8/3/2023 JZ, adapted Natalia's script to batch process Corelets tracking data
% need to add pairwiseMSD function to path

clear
clc
%close all

% Get file directions:
directory = 'E:\Dropbox\Dropbox (Princeton)\Princeton Brangwynne Lab\Data\2022-2023\20230531_GEM_squeezing\spots\';
addpath(directory)
addpath('E:\Dropbox\Dropbox (Princeton)\Princeton Brangwynne Lab\Data\matlab\msd\corelet pairwise\')
cd(directory)
dirspots = dir([directory '*back.csv']);
fnames = {};
numfiles = length(dirspots);
[fnames{1:numfiles}] = dirspots(:).name;

figure

c=0;
for fnum = 1:numfiles
    output = pairwiseMSD(fnames{fnum}); 
    for num = 1:length(output)
        c=c+1;
        table_tmp = output{num};
        msd_all(:,c) = table_tmp.MSD; 
        loglog(table_tmp.Lagtime,table_tmp.MSD,'LineWidth',0.75,'Color',[0.85,0.85,0.85, 0.2])
        hold on
    end
end

loglog(table_tmp.Lagtime,mean(msd_all,2),'LineWidth',3,'Color',[80/256, 148/256, 250/256, 0.3])
T_averaged = table(table_tmp.Lagtime, mean(msd_all,2),'VariableNames', {'lagtime','msd'});
writetable(T_averaged, sprintf('DZNep_PairwiseMSD_confined_back_avg.csv'));

xlim([0.06 1])
ylabel('Pairwise MSD (\mum^{2})')
xlabel('\tau (s)')
set(gca,'linewidth',1)
ax = gca; 
ax.FontSize = 20;
set(gca,'TickLength',[0.02, 0.02])
set(gcf,'Position',[100 100 400 400])

%fit diffusive exponent and diffusion coefficient
tdata_full = T_averaged.lagtime;
tdata = tdata_full(1:20);
ndata_full = T_averaged.msd;
ndata = ndata_full(1:20);

p = fit(tdata, ndata, 'power1');
Dapp_avg=p.a/4;
alpha_avg=p.b;

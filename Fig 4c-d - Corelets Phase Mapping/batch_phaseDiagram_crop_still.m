%12/16/23 Updated version for batch process for plotting phase diagram
%JZ
clear
clc
close

tic;

dir_n = '/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/2022-2023/20231210 1211 1221 1229 DZNep Corelet Activation/CropsNoPS/';
dir_y = '/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/2022-2023/20231210 1211 1221 1229 DZNep Corelet Activation/CropsYesPS/';
% dir_y='/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/2022-2023/FCS Dynamic Profiler/Corelet Activation Tifs/curr/single crops/Yes/';
% dir_n='/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/2022-2023/FCS Dynamic Profiler/Corelet Activation Tifs/curr/single crops/No/';
% dir_y = '/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/2022-2023/20231021 all corelets overnights optimization/20231104 corelet overnight yPS cells 15-sec activation/pre-overnight tests/RearEndCoreletsExamples/MatlabInput/Cell2/';
% dir_n = '/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/2022-2023/20231021 all corelets overnights optimization/20231104 corelet overnight yPS cells 15-sec activation/pre-overnight tests/RearEndCoreletsExamples/MatlabInput/Cell2/';

addpath(dir_y)
addpath(dir_n)

dirobj_y=dir([dir_y '*.tif']);
dirobj_n=dir([dir_n '*.tif']);

fnames_y={};
fnames_n={};
[fnames_y{1:length(dirobj_y),1}]=dirobj_y(:).name; %create cell array containing file names
[fnames_n{1:length(dirobj_n),1}]=dirobj_n(:).name; %create cell array containing file names
numfiles=length(dirobj_y)+length(dirobj_n);

full_fnames_y = strcat(dir_y,fnames_y);
full_fnames_n = strcat(dir_n,fnames_n);

full_fnames_all = [full_fnames_y; full_fnames_n];

concatInverseValence = [];
concatCoreConc = [];
concatPS = [];
 
for fnum = 1:numfiles
    preact_filename = full_fnames_all{fnum};
    BG_GFP = 100;
    BG_mCherry = 100;
    [inverse_valence_temp, core_conc_temp, call_PS_temp] = callPhaseSepStill(preact_filename, BG_GFP, BG_mCherry);
    inverse_valence{fnum,1}=inverse_valence_temp; 
    core_conc{fnum,1}=core_conc_temp;
    call_PS{fnum,1}=call_PS_temp;
    concatInverseValence = [concatInverseValence; inverse_valence{fnum,1}];
    concatCoreConc = [concatCoreConc; core_conc{fnum,1}];
    concatPS = [concatPS; call_PS{fnum,1}];
end

%% Get phase boundary by SVM classification and plot!
addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/colorspace/colorspace/')
addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/cbrewer2/')

log2val = log2(concatInverseValence);

figure
scatter(concatCoreConc(concatPS==1),log2val(concatPS==1),'MarkerFaceColor','b','MarkerFaceAlpha',0.4,'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha',0,'SizeData',25)
hold on
scatter(concatCoreConc(concatPS==0),log2val(concatPS==0),'MarkerFaceColor','m','MarkerFaceAlpha',0.2,'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha',0,'SizeData',25)
hold on
xlabel('[core] (\muM)')
ylabel('Valence^{-1}')
xlim([0 1])
y_values = [-7, -6, -5, -4, -3, -2, -1, 0];
set(gca, 'YTick', y_values, 'YTickLabel', {'1/128', '1/64','1/32', '1/16', '1/8', '1/4', '1/2', '1'});
set(gca, 'YMinorTick', 'off')
set(gca, 'LineWidth', 1.5)
ylim([-7 0])
set(gca,'FontSize',30)
set(gcf, 'Position',[100 100 400 400])
axis square

x = [concatCoreConc, log2val];
y = concatPS;
svmModel = fitcsvm(x, y, 'KernelFunction', 'polynomial', 'PolynomialOrder', 2);

% Generate a grid of values over your data range
[xGrid, yGrid] = meshgrid(linspace(0,1,300), linspace(-7,0,300));
XGrid = [xGrid(:), yGrid(:)];

% Predict labels for the grid
predictions = predict(svmModel, XGrid);

% Reshape the predictions to the grid shape
decisionMap = reshape(predictions, size(xGrid));

% Plot the decision boundary
contourf(xGrid, yGrid, decisionMap, 1, '--', 'LineWidth', 2.5, 'Color',[0.75 0.75 0.75],  'FaceAlpha',0.25);
color_index = cbrewer2('Pastel1',2);
colormap(color_index)

% %% Comment this out if not wish to plot single points
% log2val = log2(concatInverseValence./1.362);
% 
% hold on
% scatter(concatCoreConc(concatPS==1)./1.362,log2val(concatPS==1),'MarkerFaceColor','b','MarkerFaceAlpha',1,'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha',0,'SizeData',150)
% hold on
% scatter(concatCoreConc(concatPS==0)./1.362,log2val(concatPS==0),'MarkerFaceColor','m','MarkerFaceAlpha',1,'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha',0,'SizeData',150)
% hold on
% y_values = [-8, -7, -6, -5, -4];
% set(gca, 'YTick', y_values, 'YTickLabel', {'1/256', '1/128', '1/64','1/32', '1/16'});
% ylim([-8 -4])
% xlim([0 0.3])

%%
toc;
disp('Elapsed time: '); disp(toc);

%% generate phase separating data points
function [inverse_valence, core_conc, call_PS] = callPhaseSepStill(preact_filename, BG_GFP, BG_mCherry)
core_pre_uint16=imread(preact_filename,2)-BG_GFP;
core_pre_double=double(core_pre_uint16);
core_pre_double(core_pre_double<0)=0;
IDR_pre_uint16=imread(preact_filename,1)-BG_mCherry;
IDR_pre_double=double(IDR_pre_uint16);
IDR_pre_double(IDR_pre_double<0)=0;

core_filt=imgaussfilt(core_pre_double);

core_threshmat=imfilter(core_pre_double,fspecial('disk',100));
nucl_threshmat=imfilter(core_pre_double,fspecial('disk',60));
mask_core=imfill(bwareafilt(imbinarize(core_filt,core_threshmat),[1e3 8e4]),'holes');

mask_nucleoli = bwareafilt(imcomplement(bwareafilt(imbinarize(core_filt,nucl_threshmat),[4e3 5e4])),[200,4000]);
mask_adjusted = mask_core.*(~mask_nucleoli);
mask_adjusted = logical(imerode(mask_adjusted,strel('disk',2)));
mask_final = bwareafilt(mask_adjusted, [1e3, 8e4]);

core_conc = [];
IDR_conc = [];
[labelmat, number_nucleus] = bwlabel(mask_final,8);

for n = 1:number_nucleus
    core_conc(n,1) = (mean(core_pre_double(labelmat==n),'all') .* 44.6/1000)/24; %convert to micro-molar
    IDR_conc(n,1) = mean(IDR_pre_double(labelmat==n),'all') .* 44.6/1000/0.633; % convert to micro-molar
end
inverse_valence = core_conc./IDR_conc; 

call_PS = zeros(number_nucleus,1);
if contains(preact_filename,'Yes')
    call_PS(:) = 1;
end

end
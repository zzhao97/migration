clear
clc
close all

% Get file directions:
directory = 'E:\Dropbox\Princeton Dropbox\Jessica Zhao\Princeton Brangwynne Lab\Data\Confined Migration Paper\Nature Comm Revision\New Experiments\H2B Heterogeneity\tiff\';
addpath(directory)
addpath('E:\Dropbox\Princeton Dropbox\Jessica Zhao\Princeton Brangwynne Lab\Data\matlab')
cd(directory)
dirspots = dir([directory '*.tif']);
fnames = {};
numfiles = length(dirspots);
[fnames{1:numfiles}] = dirspots(:).name;

H2B_cov_array = [];
names_array = [];

for fnum = 1:numfiles
    filename = fnames{fnum};
    [cov_temp, number_nucleus_temp] = getCoV(filename);
    names_arr_temp = repmat(filename,number_nucleus_temp,1);
    H2B_cov_array = [H2B_cov_array; cov_temp];
    names_array = [names_array; cellstr(names_arr_temp)];
end

% save output table as csv file
csv_table = table(names_array, H2B_cov_array, 'VariableNames',{'filename','H2B_Heterogeneity'});
writetable(csv_table,'H2B_Heterogeneity_all_20240818.csv');


function [cov, number_nucleus] = getCoV(filename)
H2B_img = imread(filename);

% background subtraction
H2B_no_bg = H2B_img - mean(H2B_img(1:10, 1:10),'all');
H2B_no_bg(H2B_no_bg<0) = 0;

H2B_double_img = double(H2B_no_bg);
H2B_filt=imgaussfilt(H2B_double_img,2);

% segment nucleus area
beta = 0.25;
se = strel('disk',3);

thresh=mean(H2B_filt,'all')+beta*std(H2B_filt,0,'all');
mask_nuc=imclearborder(bwareafilt(imfill(imbinarize(H2B_filt,thresh),'holes'),[1e4 5e4]));

[labelmat, number_nucleus] = bwlabel(mask_nuc,8);
%figure; imagesc(labelmat); axis equal;

for n = 1:number_nucleus
    % calculate coefficient of variation as metric for H2B heterogeneity
    all_H2B_pixels = H2B_filt(labelmat==n);
    cov(n,1) = std(all_H2B_pixels,0,'all')./mean(all_H2B_pixels,'all');
end

end
% Confined migration paper JZ 2024
%
%
% - Used in analyzing perinucleolar heterochromatin intensity
% - Reads multi-color tif stacks with H2B and nucleoli channels
% - Input tif files are manually cropped so only contain trailing halves
% - H2B channel used to segment trailing half nucleus
% then calculate size of trailing half for every time point
% - Inputs tifs are front three different independent replicates
%
%
clear; 
close all; 
clc;

output_dir = '/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/2024/20240630_NPM1_perinucleolar_chromatin/backAreaOutputs/';
cd(output_dir)

directory_data='/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/2024/20240630_NPM1_perinucleolar_chromatin/back tif/';

addpath(directory_data)

dir_data = dir([directory_data '*ctrl*.tif']);
fnames_data={};
[fnames_data{1:length(dir_data)}]=dir_data(:).name;

for fnum = 1:length(dir_data) % change folder if looking at another dataset
    filename = fnames_data{fnum}; % change folder if looking at another dataset
    back_area_temp = getBackArea(filename);
    Tdata = table(repmat(filename,[length(back_area_temp) 1]), back_area_temp, 'VariableNames',{'Filename','BackArea'});
    writetable(Tdata, sprintf('BackArea_%s.csv', fnames_data{fnum})) % change filename and folder
end

% perinucleolar heterochromatin condensation
function [back_area] = getBackArea(filename)

% read file and background subtraction
size_img=size(imfinfo(filename),1);
ncolors=3;
nframes=size_img/ncolors;

H2B_img=zeros([size(imread(filename),1)+4 size(imread(filename),2) nframes]);

for k=1:nframes
    H2B_img(5:end,:,k)=imread(filename,3*k-2);
end

background_H2B=2;
H2B_img=H2B_img-background_H2B;

% segment nucleus
for k=1:nframes
    H2B_img_filt(:,:,k)=imgaussfilt(H2B_img(:,:,k),0.25);
    thresh_H2B=mean(H2B_img_filt(:,:,k),'all')+0.25*std(H2B_img_filt(:,:,k),0,'all');
    mask_H2B(:,:,k)=imclearborder(bwareafilt(imfill(imbinarize(H2B_img_filt(:,:,k),thresh_H2B),'holes'),[2e3 9e4]));
    back_area(k,1)=sum(mask_H2B(:,:,k),[1 2]).*0.085^2; %um^2
end

end
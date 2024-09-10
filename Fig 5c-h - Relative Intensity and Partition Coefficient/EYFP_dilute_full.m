% Confined Migration Paper 2024
%
%
% This code is for quantifying BACK HALF of the squeezed nucleus
% to quantify: 
% EYFP_CoV: coefficient of variation - how heterogeneous the signal dist is 
% EYFP_mat_b: material in the back - total signal in EYFP channel  
% area_b: area of the back half
% update on 3/25/24

clear; 
clc; 
close all;

directory = '/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/2022-2023/20231210_EYFP_H2B-miRFP_overnight/';

addpath(directory)
addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/MDA MB 231 confined migration features/EYFP rear-whole distribution scripts/')
cd(directory)
dirobj=dir([directory '*-f.tif']);%find all tifs in the directory
[fnames{1:length(dirobj)}]=dirobj(:).name; %create cell array containing file names
numfiles=length(dirobj);
area_b_array = [];
filenamelist = {}; collapse = [];

for f = 1:numfiles
    filename = fnames{f};
    if contains(filename, 'only')
        [area_b_temp, nframes] = getFullMeasurementsSingle(filename);
    else
        [area_b_temp, nframes] = getFullMeasurementsDual(filename);
    end
    EYFP_area_b{f,1}=area_b_temp;
    area_b_array = [area_b_array; EYFP_area_b{f,1}];
    filenamelist{f,1}=repmat(cellstr(filename),[nframes,1]);
    collapse = [collapse; string(filenamelist{f,1})];
end

Tdata=table(collapse,area_b_array,'VariableNames',{'Filename','AreaFull'}) %generate output table
writetable(Tdata,[directory 'SRRM1_full_output_partition.csv'])


%% EYFP only cells (lacks H2B-miRFP as nucleus marker)
function [area_b, nframes] = getFullMeasurementsSingle(filename)
size_img=size(imfinfo(filename),1);
ncolors=3;
nframes=size_img/ncolors;


%pre-allocate for speed
EYFP_img_orig=zeros([size(imread(filename),1) size(imread(filename),2) nframes]);
EYFP_img=EYFP_img_orig;


%read in multistack
for k=1:nframes
    EYFP_img_orig(:,:,k)=imread(filename,3*k-2);
    EYFP_img=double(EYFP_img_orig);
end


%subtract background for both H2B and protein images
EYFP_img=EYFP_img-15;
EYFP_img(EYFP_img<0)=0;


%segment nucleus area
se = strel('disk',3);
beta = 0;

for k=1:nframes
    EYFP_filt(:,:,k)=imgaussfilt(EYFP_img(:,:,k),1.5);
    thresh=mean(EYFP_filt(:,:,k),'all')+beta*std(EYFP_filt(:,:,k),0,'all');
    mask_EYFP(:,:,k)=bwareafilt(imbinarize(imgaussfilt(EYFP_filt(:,:,k),2),thresh),[1e3 5e4]);
    mask_nuc(:,:,k)=imfill(imclearborder(mask_EYFP(:,:,k)),'holes');
%     figure
%     imagesc(mask_nuc(:,:,k)), axis equal
    area_b(k,1)=sum(mask_nuc(:,:,k),[1 2]).*0.085^2; %um^2
end
end

%% H2B-miRFP expressing cells, using it as a nucleus marker for segmentation
function [area_b, nframes] = getFullMeasurementsDual(filename)
size_img=size(imfinfo(filename),1);
ncolors=3;
nframes=size_img/ncolors;


%pre-allocate for speed
EYFP_img_orig=zeros([size(imread(filename),1) size(imread(filename),2) nframes]);
H2B_img_orig=EYFP_img_orig;

EYFP_img=EYFP_img_orig;
H2B_img=EYFP_img_orig;


%read in multistack
for k=1:nframes
    EYFP_img_orig(:,:,k)=imread(filename,3*k-2);
    H2B_img_orig(:,:,k)=imread(filename,3*k-1);
    H2B_img=double(H2B_img_orig);
    EYFP_img=double(EYFP_img_orig);
end


%subtract background for both H2B and protein images
H2B_img=H2B_img-0;
EYFP_img=EYFP_img-15;
H2B_img(H2B_img<0)=0;
EYFP_img(EYFP_img<0)=0;


%segment nucleus area
se = strel('disk',3);
beta = 0.5;

for k=1:nframes
    H2B_filt(:,:,k)=imgaussfilt(H2B_img(:,:,k),2);
    EYFP_filt(:,:,k)=imgaussfilt(EYFP_img(:,:,k),1.5);
    thresh=mean(H2B_filt(:,:,k),'all')+beta*std(H2B_filt(:,:,k),0,'all');
    mask_H2B(:,:,k)=imclearborder(bwareafilt(imbinarize(H2B_filt(:,:,k),thresh),[1e3 8e4]));
    mask_nuc(:,:,k)=imclearborder(bwareafilt(imfill(imbinarize(H2B_filt(:,:,k),thresh),'holes'),[1e3 8e4]));
 
    area_b(k,1)=sum(mask_nuc(:,:,k),[1 2]).*0.085^2; %um^2

end

end
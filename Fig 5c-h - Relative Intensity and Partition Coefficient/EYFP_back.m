% Confined Migration Paper 2024
%
%
% This code is for quantifying BACK HALF of the squeezed nucleus
% to quantify: 
% EYFP_CoV: coefficient of variation - how heterogeneous the signal dist is 
% EYFP_mat_b: material in the back - total signal in EYFP channel  
% area_b: area of the back half
% 
% 

clear; 
clc; 
close all;

directory='/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/2022-2023/20231210_EYFP_H2B-miRFP_overnight/';

addpath(directory)
addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/MDA MB 231 confined migration features/EYFP rear-whole distribution scripts/')
cd(directory)
dirobj=dir([directory '*-b.tif']);%find all tifs in the directory
fnames={}; 
[fnames{1:length(dirobj)}]=dirobj(:).name; %create cell array containing file names
numfiles=length(dirobj);
EYFP_cov_array = [];
EYFP_mat_b_array= []; 
area_b_array = [];
filenamelist = {}; collapse = [];

for f = 1:numfiles
    filename = fnames{f};
    if contains(filename, 'only')
        [cov_temp, mat_b_temp, area_b_temp, nframes] = getBackMeasurementsSingle(filename);
    else
        [cov_temp, mat_b_temp, area_b_temp, nframes] = getBackMeasurementsDual(filename);
    end
    EYFP_cov{f,1}=cov_temp; 
    EYFP_mat_b{f,1}=mat_b_temp;
    EYFP_area_b{f,1}=area_b_temp;
    EYFP_cov_array = [EYFP_cov_array; EYFP_cov{f,1}];
    EYFP_mat_b_array = [EYFP_mat_b_array; EYFP_mat_b{f,1}];
    area_b_array = [area_b_array; EYFP_area_b{f,1}];
    filenamelist{f,1}=repmat(cellstr(filename),[nframes,1]);
    collapse = [collapse; string(filenamelist{f,1})];
end

Tdata=table(collapse, EYFP_cov_array,EYFP_mat_b_array,area_b_array,'VariableNames',{'Filename','CoV','MaterialBack','AreaBack'}) %generate output table
writetable(Tdata,[directory 'EYFP_only_back_output_data.csv'])


%% EYFP-channel only cells (lacks H2B-miRFP as nucleus marker)
function [cov, mat_b, area_b, nframes] = getBackMeasurementsSingle(filename)
size_img=size(imfinfo(filename),1);
ncolors=3;
nframes=size_img/ncolors;


%pre-allocate for speed
EYFP_img_orig=zeros([size(imread(filename),1)+4 size(imread(filename),2) nframes]);
%this will initialize a matrix with first dimension 4 pixels greater so when
%imclearborder runs it doesn't get rid of the half nucleus
EYFP_img=EYFP_img_orig;


%read in multistack
for k=1:nframes
    EYFP_img_orig(5:end,:,k)=imread(filename,3*k-2);
    %get rid of dark border in Fiji cropping. not sure why this occurs
    EYFP_img=double(EYFP_img_orig);
end


%subtract background for both H2B and protein images
EYFP_img=EYFP_img-0;
EYFP_img(EYFP_img<0)=0;


%segment nucleus area
se = strel('disk',3);
beta = 0;

for k=1:nframes
    EYFP_filt(:,:,k)=imgaussfilt(EYFP_img(:,:,k),1.5);
    thresh=mean(EYFP_filt(:,:,k),'all')+beta*std(EYFP_filt(:,:,k),0,'all');
    mask_EYFP(:,:,k)=bwareafilt(imbinarize(imgaussfilt(EYFP_filt(:,:,k),2),thresh),[8e2 5e4]);
    mask_nuc(:,:,k)=imfill(imclearborder(mask_EYFP(:,:,k)),'holes');
%     figure
%     imagesc(mask_nuc(:,:,k)), axis equal
    area_b(k,1)=sum(mask_nuc(:,:,k),[1 2]).*0.085^2; %um^2
    
    mat_b(k,1)=sum(EYFP_filt(:,:,k).*mask_nuc(:,:,k),'all');
%     disp(RNPS1_mat_b(k,1));
     
    nucl_threshmat(:,:,k)=imfilter(EYFP_img(:,:,k),fspecial('disk',60));
    mask_nucleoli(:,:,k) = bwareafilt(imcomplement(bwareafilt(imbinarize(EYFP_filt(:,:,k),nucl_threshmat(:,:,k)),[8e2 5e4])),[200,3800]);
    mask_adjusted(:,:,k) = mask_nuc(:,:,k).*(~mask_nucleoli(:,:,k));
    mask_adjusted(:,:,k) = logical(imclearborder(imerode(mask_adjusted(:,:,k),se)));
    local_mask(:,:,k) = bwareafilt(logical(mask_adjusted(:,:,k)), [8e2, 5e4]);
%     figure
%     imagesc(local_mask(:,:,k)), axis equal
    
    cov(k,1)=std(EYFP_filt(:,:,k).*local_mask(:,:,k),0,'all')./mean(EYFP_filt(:,:,k).*local_mask(:,:,k),'all');

end
end

%% H2B-miRFP expressing cells, using it as a nucleus marker for segmentation
function [cov, mat_b, area_b, nframes] = getBackMeasurementsDual(filename)
size_img=size(imfinfo(filename),1);
ncolors=3;
nframes=size_img/ncolors;


%pre-allocate for speed
EYFP_img_orig=zeros([size(imread(filename),1)+4 size(imread(filename),2) nframes]);
H2B_img_orig=EYFP_img_orig;
%this will initialize a matrix with first dimension 4 pixels greater so when
%imclearborder runs it doesn't get rid of the half nucleus
EYFP_img=EYFP_img_orig;
H2B_img=EYFP_img_orig;


%read in multistack
for k=1:nframes
    EYFP_img_orig(5:end,:,k)=imread(filename,3*k-2);
    H2B_img_orig(5:end,:,k)=imread(filename,3*k-1);
    H2B_img=double(H2B_img_orig);
    EYFP_img=double(EYFP_img_orig);
end


%subtract background for both H2B and protein images
H2B_img=H2B_img-0;
EYFP_img=EYFP_img-0;
H2B_img(H2B_img<0)=0;
EYFP_img(EYFP_img<0)=0;


%segment nucleus area
se = strel('disk',3);
beta = 0.5;

for k=1:nframes
    H2B_filt(:,:,k)=imgaussfilt(H2B_img(:,:,k),2);
    EYFP_filt(:,:,k)=imgaussfilt(EYFP_img(:,:,k),1.5);
    thresh=mean(H2B_filt(:,:,k),'all')+beta*std(H2B_filt(:,:,k),0,'all');
    mask_H2B(:,:,k)=imclearborder(bwareafilt(imbinarize(H2B_filt(:,:,k),thresh),[8e2 5e4]));
    mask_nuc(:,:,k)=imclearborder(bwareafilt(imfill(imbinarize(H2B_filt(:,:,k),thresh),'holes'),[8e2 5e4]));
%     figure
%     imagesc(mask_nuc(:,:,k)), axis equal
    area_b(k,1)=sum(mask_nuc(:,:,k),[1 2]).*0.085^2; %um^2
    
    mat_b(k,1)=sum(EYFP_filt(:,:,k).*mask_nuc(:,:,k),'all');
%     disp(RNPS1_mat_b(k,1));
     
    nucl_threshmat(:,:,k)=imfilter(EYFP_img(:,:,k),fspecial('disk',60));
    mask_nucleoli(:,:,k) = bwareafilt(imcomplement(bwareafilt(imbinarize(EYFP_filt(:,:,k),nucl_threshmat(:,:,k)),[8e2 5e4])),[200,3800]);
    mask_adjusted(:,:,k) = mask_nuc(:,:,k).*(~mask_nucleoli(:,:,k));
    mask_adjusted(:,:,k) = logical(imclearborder(imerode(mask_adjusted(:,:,k),se)));
    local_mask(:,:,k) = bwareafilt(logical(mask_adjusted(:,:,k)), [1e3, 5e4]);
%     figure
%     imagesc(mask_final(:,:,k)), axis equal
    
    cov(k,1)=std(EYFP_filt(:,:,k).*local_mask(:,:,k),0,'all')./mean(EYFP_filt(:,:,k).*local_mask(:,:,k),'all');

end

end
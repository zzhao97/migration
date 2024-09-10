% This code is for quantifying BACK HALF of the squeezed nucleus
% to quantify: 
% EYFP_CoV: coefficient of variation - how heterogeneous the signal dist is 
% EYFP_mat_b: material in the back - total signal in EYFP channel  
% area_b: area of the back half
% update on 11/4/23 

clear; 
clc; 
close all;

directory='/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/2022-2023/20231227_SRRM1_H2B_migration_duplicates_all/';
%directory='/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/20230425_speckle_dynamics_migration/';
%directory='/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/20231016 all EYFP constructs migration/20231031_DAXX-EYFP_high_crops/';
%directory='/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/20231016 all EYFP constructs migration/20231018_SART1-EYFP_H2B-miRFP_overnight/';
%directory='/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/20231016 all EYFP constructs migration/20231016_RNPS1-EYFP_H2B-miRFP_highExp_8uL_overnight/';

addpath(directory)
addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/MDA MB 231 confined migration features/H2B rear-front scripts/')
cd(directory)
dirobj=dir([directory '*-b.tif']);%find all tifs in the directory
fnames={}; 
[fnames{1:length(dirobj)}]=dirobj(:).name; %create cell array containing file names
numfiles=length(dirobj);
H2B_cov_array = [];
H2B_mat_b_array= []; 
area_b_array = [];
filenamelist = {}; collapse = [];

for f = 1:numfiles
    filename = fnames{f};
    if ~contains(filename,'only') % skip files with only protein channel
        [cov_temp, mat_b_temp, area_b_temp, nframes] = getBackMeasurementsDual(filename);
        H2B_cov{f,1}=cov_temp./cov_temp(1);
        H2B_mat_b{f,1}=mat_b_temp;
        H2B_area_b{f,1}=area_b_temp;
        H2B_cov_array = [H2B_cov_array; H2B_cov{f,1}];
        H2B_mat_b_array = [H2B_mat_b_array; H2B_mat_b{f,1}];
        area_b_array = [area_b_array; H2B_area_b{f,1}];
        filenamelist{f,1}=repmat(cellstr(filename),[nframes,1]);
        collapse = [collapse; string(filenamelist{f,1})];
    end
end

Tdata=table(collapse, H2B_cov_array,H2B_mat_b_array,area_b_array,'VariableNames',{'Filename','CoV','MaterialBack','AreaBack'}) %generate output table
writetable(Tdata,[directory 'H2B_back_output_data.csv'])

%% 
function [cov, mat_b, area_b, nframes] = getBackMeasurementsDual(filename)
size_img=size(imfinfo(filename),1);
ncolors=3;
nframes=size_img/ncolors;

%pre-allocate for speed
H2B_img_orig=zeros([size(imread(filename),1)+4 size(imread(filename),2) nframes]);
%this will initialize a matrix with first dimension 4 pixels greater so when
%imclearborder runs it doesn't get rid of the half nucleus
H2B_img=H2B_img_orig;

%read in multistack
for k=1:nframes
    H2B_img_orig(5:end,:,k)=imread(filename,3*k-1);
    H2B_img=double(H2B_img_orig);
end

%subtract background for both H2B and protein images
H2B_img=H2B_img-4;
H2B_img(H2B_img<0)=0;

%segment nucleus area
se = strel('disk',3);
beta = 0.5;

for k=1:nframes
    H2B_filt(:,:,k)=imgaussfilt(H2B_img(:,:,k),2);
    thresh=mean(H2B_filt(:,:,k),'all')+beta*std(H2B_filt(:,:,k),0,'all');
    mask_H2B(:,:,k)=imclearborder(bwareafilt(imbinarize(H2B_filt(:,:,k),thresh),[8e2 5e4]));
    mask_nuc(:,:,k)=imclearborder(bwareafilt(imfill(imbinarize(H2B_filt(:,:,k),thresh),'holes'),[8e2 5e4]));

    area_b(k,1)=sum(mask_nuc(:,:,k),[1 2]).*0.085^2; %um^2
    
    mat_b(k,1)=sum(H2B_filt(:,:,k).*mask_nuc(:,:,k),'all');
    
    cov(k,1)=std(H2B_filt(:,:,k).*mask_nuc(:,:,k),0,'all')./mean(H2B_filt(:,:,k).*mask_nuc(:,:,k),'all');

end

end
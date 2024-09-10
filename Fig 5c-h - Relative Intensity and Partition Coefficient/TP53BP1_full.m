% This code is for quantifying BACK HALF of the squeezed nucleus
% to quantify: 
% TP53BP1_CoV: coefficient of variation - how heterogeneous the signal dist is 
% TP53BP1_mat_b: material in the back - total signal in TP53BP1 channel  
% area_b: area of the back half
% update on 11/4/23 

clear; 
clc; 
close all;

directory='/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/20230111 53BP1 concentration/no squeeze/';

addpath(directory)
addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/MDA MB 231 confined migration features/53BP1 rear-whole scripts/')

cd(directory)
dirobj=dir([directory '*-whole.tif']);%find all tifs in the directory
fnames={}; 
[fnames{1:length(dirobj)}]=dirobj(:).name; %create cell array containing file names
numfiles=length(dirobj);
TP53BP1_cov_array = [];
TP53BP1_mat_b_array= []; 
area_b_array = [];
filenamelist = {}; collapse = [];

for f = 1:numfiles
    filename = fnames{f};
    [cov_temp, mat_b_temp, area_b_temp, nframes] = getFullMeasurementsDual(filename);
    TP53BP1_cov{f,1}=cov_temp; 
    TP53BP1_mat_b{f,1}=mat_b_temp;
    TP53BP1_area_b{f,1}=area_b_temp;
    TP53BP1_cov_array = [TP53BP1_cov_array; TP53BP1_cov{f,1}];
    TP53BP1_mat_b_array = [TP53BP1_mat_b_array; TP53BP1_mat_b{f,1}];
    area_b_array = [area_b_array; TP53BP1_area_b{f,1}];
    filenamelist{f,1}=repmat(cellstr(filename),[nframes,1]);
    collapse = [collapse; string(filenamelist{f,1})];
end

Tdata=table(collapse, TP53BP1_cov_array,TP53BP1_mat_b_array,area_b_array,'VariableNames',{'Filename','CoV','MaterialFull','AreaFull'}) %generate output table
writetable(Tdata,[directory '20231123_no_squeeze_whole_output_data.csv'])

%% H2B-GFP expressing cells, using it as a nucleus marker for segmentation
function [cov, mat_b, area_b, nframes] = getFullMeasurementsDual(filename)
size_img=size(imfinfo(filename),1);
ncolors=3;
nframes=size_img/ncolors;


%pre-allocate for speed
TP53BP1_img_orig=zeros([size(imread(filename),1) size(imread(filename),2) nframes]);
H2B_img_orig=TP53BP1_img_orig;

TP53BP1_img=TP53BP1_img_orig;
H2B_img=TP53BP1_img_orig;


%read in multistack
for k=1:nframes
    TP53BP1_img_orig(:,:,k)=imread(filename,3*k-2);
    H2B_img_orig(:,:,k)=imread(filename,3*k-1);
    H2B_img=double(H2B_img_orig);
    TP53BP1_img=double(TP53BP1_img_orig);
end


%subtract background for both H2B and protein images
H2B_img=H2B_img-4;
TP53BP1_img=TP53BP1_img-15;
H2B_img(H2B_img<0)=0;
TP53BP1_img(TP53BP1_img<0)=0;


%segment nucleus area
se = strel('disk',3);
beta = 0.5;

for k=1:nframes
    H2B_filt(:,:,k)=imgaussfilt(H2B_img(:,:,k),2);
    TP53BP1_filt(:,:,k)=imgaussfilt(TP53BP1_img(:,:,k),1.5);
    thresh=mean(H2B_filt(:,:,k),'all')+beta*std(H2B_filt(:,:,k),0,'all');
    mask_H2B(:,:,k)=imclearborder(bwareafilt(imbinarize(H2B_filt(:,:,k),thresh),[1e3 5e4]));
    mask_nuc(:,:,k)=imclearborder(bwareafilt(imfill(imbinarize(H2B_filt(:,:,k),thresh),'holes'),[1e3 5e4]));
%     figure
%     imagesc(mask_nuc(:,:,k)), axis equal
    area_b(k,1)=sum(mask_nuc(:,:,k),[1 2]).*0.085^2; %um^2
    
    mat_b(k,1)=sum(TP53BP1_filt(:,:,k).*mask_nuc(:,:,k),'all');
%     disp(TP53BP1_mat_b(k,1));
     
    nucl_threshmat(:,:,k)=imfilter(TP53BP1_img(:,:,k),fspecial('disk',60));
    mask_nucleoli(:,:,k) = bwareafilt(imcomplement(bwareafilt(imbinarize(TP53BP1_filt(:,:,k),nucl_threshmat(:,:,k)),[1e3 5e4])),[200,3800]);
    mask_adjusted(:,:,k) = mask_nuc(:,:,k).*(~mask_nucleoli(:,:,k));
    mask_adjusted(:,:,k) = logical(imclearborder(imerode(mask_adjusted(:,:,k),se)));
    local_mask(:,:,k) = bwareafilt(logical(mask_adjusted(:,:,k)), [1e3, 5e4]);
%     figure
%     imagesc(mask_final(:,:,k)), axis equal
    
    cov(k,1)=std(TP53BP1_filt(:,:,k).*local_mask(:,:,k),0,'all')./mean(TP53BP1_filt(:,:,k).*local_mask(:,:,k),'all');

end

end
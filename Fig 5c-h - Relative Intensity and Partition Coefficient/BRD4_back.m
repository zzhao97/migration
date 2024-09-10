% This code is for quantifying back half of squeezed nucleus
% Goal: use BRD4 channel to segment nucleus for getting BRD4
% signal only in the nucleoplasm
% mat_b: BRD4 material in the back - total signal of BRD4 in nucleoplasm
% area_b: area of the half nucleus
%
% update on 3/6/24 (BRD4 channel only)

clear;
clc;
close all;

input_directory = '/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/2024/20240228_BRD4FL-miRFP_migration_pilot/single crops/';
addpath(input_directory)
cd(input_directory)

addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/MDA MB 231 confined migration features/EYFP rear-whole distribution scripts/')

dirobj=dir([input_directory '*back.tif']);%find all tifs in the directory
fnames={};
[fnames{1:length(dirobj)}]=dirobj(:).name; %create cell array containing file names
numfiles=length(dirobj);
BRD4_cov_array = [];
BRD4_mat_b_array= [];
area_b_array = [];
filenamelist = {}; collapse = [];

for f = 1:numfiles
    filename = fnames{f};
    [mat_b_temp, area_b_temp, nframes] = getBackMeasurementsSingle(filename);
    BRD4_mat_b{f,1}=mat_b_temp;
    BRD4_area_b{f,1}=area_b_temp;
    BRD4_mat_b_array = [BRD4_mat_b_array; BRD4_mat_b{f,1}];
    area_b_array = [area_b_array; BRD4_area_b{f,1}];
    filenamelist{f,1}=repmat(cellstr(filename),[nframes,1]);
    collapse = [collapse; string(filenamelist{f,1})];
end

Tdata=table(collapse, BRD4_mat_b_array, area_b_array,'VariableNames',{'Filename','MaterialFull','AreaFull'}) %generate output table
writetable(Tdata,[input_directory 'BRD4_back_output_data_20240306.csv'])


%% BRD4-channel only cells
function [mat_b, area_b, nframes] = getBackMeasurementsSingle(filename)
size_img=size(imfinfo(filename),1);
ncolors=2;
nframes=size_img/ncolors;


%pre-allocate for speed
BRD4_img_orig=zeros([size(imread(filename),1)+4 size(imread(filename),2) nframes]);
%this will initialize a matrix with first dimension 4 pixels greater so when
%imclearborder runs it doesn't get rid of the half nucleus
BRD4_img=BRD4_img_orig;


%read in multistack
for k=1:nframes
    BRD4_img_orig(5:end,:,k)=imread(filename,2*k-1);
    BRD4_img=double(BRD4_img_orig);
end


%subtract background for BRD4 images
BRD4_img=BRD4_img-4;
BRD4_img(BRD4_img<0)=0;


%segment nucleus area
beta_dil = 0.5;
se = strel('disk',3);

for k=1:nframes
    BRD4_filt(:,:,k)=imgaussfilt(BRD4_img(:,:,k),2);
    thresh=mean(BRD4_filt(:,:,k),'all')+beta_dil*std(BRD4_filt(:,:,k),0,'all');
    mask_nuc(:,:,k)=imclearborder(bwareafilt(imfill(imbinarize(BRD4_filt(:,:,k),thresh),'holes'),[8e2 4e4]));

    nucl_threshmat(:,:,k)=imfilter(BRD4_filt(:,:,k),fspecial('disk',60));
    mask_nucleoli(:,:,k) = bwareafilt(imcomplement(bwareafilt(imbinarize(BRD4_filt(:,:,k),nucl_threshmat(:,:,k)),[8e2 5e4])),[100,3800]);
    mask_adjusted(:,:,k) = mask_nuc(:,:,k).*(~mask_nucleoli(:,:,k));
    mask_adjusted(:,:,k) = logical(imclearborder(imerode(mask_adjusted(:,:,k),se)));
    local_mask(:,:,k) = bwareafilt(logical(mask_adjusted(:,:,k)), [8e2, 5e4]);

    area_b(k,1)=sum(mask_nuc(:,:,k),[1 2]).*0.085^2; %um^2

    mat_b(k,1)=sum(BRD4_filt(:,:,k).*local_mask(:,:,k),'all');

    % save final mask at each frame as png, can comment this out
    % figure
    % colormap("parula")
    % imagesc(BRD4_filt(:,:,k))
    % axis equal, axis off
    % hold on
    % boundary = bwboundaries(local_mask(:,:,k));
    % for b = 1:length(boundary)
    %     plot(boundary{b}(:,2), boundary{b}(:,1),'Color','m','LineWidth',1);
    %     hold on
    % end
    % Image = getframe(gcf);
    % imwrite(Image.cdata, sprintf('BRD4_mask_%s_frame%d.png', filename, k));
    % close
end
end
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

directory = '/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/2022-2023/20231227_SRRM1_H2B_migration_duplicates_all/';

addpath(directory)
addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/MDA MB 231 confined migration features/EYFP rear-whole distribution scripts/')
cd(directory)
dirobj=dir([directory '*b.tif']);%find all tifs in the directory
fnames={};
[fnames{1:length(dirobj)}]=dirobj(:).name; %create cell array containing file names
numfiles=length(dirobj);
EYFP_PC_array = [];
EYFP_mat_b_array= [];
area_b_array = [];
filenamelist = {}; collapse = [];

for f = 1:numfiles
    filename = fnames{f};
    if contains(filename, 'only')
        [PC_temp, mat_b_temp, area_b_temp, nframes] = getBackMeasurementsSingle(filename);
    else
        [PC_temp, mat_b_temp, area_b_temp, nframes] = getBackMeasurementsDual(filename);
    end
    EYFP_PC{f,1}=PC_temp;
    EYFP_mat_b{f,1}=mat_b_temp;
    EYFP_area_b{f,1}=area_b_temp;
    EYFP_PC_array = [EYFP_PC_array; EYFP_PC{f,1}];
    EYFP_mat_b_array = [EYFP_mat_b_array; EYFP_mat_b{f,1}];
    area_b_array = [area_b_array; EYFP_area_b{f,1}];
    filenamelist{f,1}=repmat(cellstr(filename),[nframes,1]);
    collapse = [collapse; string(filenamelist{f,1})];
end

Tdata=table(collapse, EYFP_PC_array,EYFP_mat_b_array,area_b_array,'VariableNames',{'Filename','PC','MaterialBack','AreaBack'}) %generate output table
writetable(Tdata,[directory 'SRRM1_back_output_partition.csv'])


%% EYFP-channel only cells (lacks H2B-miRFP as nucleus marker)
function [PC, mat_b, area_b, nframes] = getBackMeasurementsSingle(filename)
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
EYFP_img=EYFP_img-15;
EYFP_img(EYFP_img<0)=0;


%segment nucleus area
se = strel('disk',3);
beta = 0;

for k=1:nframes
    EYFP_filt(:,:,k)=imgaussfilt(EYFP_img(:,:,k),1.5);
    thresh=mean(EYFP_filt(:,:,k),'all')+beta*std(EYFP_filt(:,:,k),0,'all');
    mask_EYFP(:,:,k)=bwareafilt(imbinarize(imgaussfilt(EYFP_filt(:,:,k),2),thresh),[8e2 5e4]);
    mask_nuc(:,:,k)=imfill(imclearborder(mask_EYFP(:,:,k)),'holes');

    area_b(k,1)=sum(mask_nuc(:,:,k),[1 2]).*0.085^2; %um^2

    mat_b(k,1)=sum(EYFP_filt(:,:,k).*mask_nuc(:,:,k),'all');

    local_mask(:,:,k) = bwareafilt(logical(mask_nuc(:,:,k)), [8e2, 5e4]);

    % segment dense phase
    circavgfilt=fspecial('disk',15);
    meanimg=imfilter(EYFP_filt(:,:,k),circavgfilt);
    imgsq=EYFP_filt(:,:,k).^2;
    meansq=imfilter(imgsq,circavgfilt);
    stdevimg=(meansq-meanimg.^2).^0.5;
    threshmat=meanimg+0.5*stdevimg;
    mask_EYFP_final(:,:,k)=bwareafilt(imfill(imbinarize(EYFP_filt(:,:,k),threshmat),'holes'),[80 800]);
    mask_EYFP_final(:,:,k)=imerode(mask_EYFP_final(:,:,k).*local_mask(:,:,k),strel('disk',1));
    mask_EYFP_final(:,:,k)=mask_EYFP_final(:,:,k).*local_mask(:,:,k);
    mask_EYFP_final(:,:,k)=bwareafilt(logical(mask_EYFP_final(:,:,k)),[80 800]);

    % segment dilute phase by dilating the dense phase binary masks
    mask_EYFP_dilate_max(:,:,k) = imdilate(mask_EYFP_final(:,:,k),strel('disk',6));
    mask_EYFP_dilate(:,:,k) = imdilate(mask_EYFP_final(:,:,k),strel('disk',2));
    dilute_mask(:,:,k) = mask_EYFP_dilate_max(:,:,k)-mask_EYFP_dilate(:,:,k);

    % save final mask at each frame as png, can comment this out
    colormap(gray)
    imagesc(EYFP_filt(:,:,k))
    axis equal, axis off
    hold on
    boundary = bwboundaries(mask_EYFP_final(:,:,k));
    for b = 1:length(boundary)
        plot(boundary{b}(:,2), boundary{b}(:,1),'Color','m','LineWidth',1);
        hold on
    end
    dilute = bwboundaries(dilute_mask(:,:,k));
    for c = 1:length(dilute)
        plot(dilute{c}(:,2), dilute{c}(:,1),'Color','b','LineWidth',1);
        hold on
    end
    Image = getframe(gcf);
    imwrite(Image.cdata, sprintf('EYFP_mask_%s_frame%d.png', filename, k));
    close

    % dense and dilute phase intensity
    EYFP_dilute_int(k,1) = mean(nonzeros(dilute_mask(:,:,k).*EYFP_filt(:,:,k)),'all');
    EYFP_dense_int(k,1) = mean(nonzeros(mask_EYFP_dilate(:,:,k).*EYFP_filt(:,:,k)),'all');

    PC(k,1) = EYFP_dense_int(k,1)./EYFP_dilute_int(k,1);

end
end

%% H2B-miRFP expressing cells, using it as a nucleus marker for segmentation
function [PC, mat_b, area_b, nframes] = getBackMeasurementsDual(filename)
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
EYFP_img=EYFP_img-15;
H2B_img(H2B_img<0)=0;
EYFP_img(EYFP_img<0)=0;


%segment nucleus area
se = strel('disk',3);
beta = 0.5;
yb = 1.5; %0 for RNPS1 and SART1

for k=1:nframes
    H2B_filt(:,:,k)=imgaussfilt(H2B_img(:,:,k),2);
    EYFP_filt(:,:,k)=imgaussfilt(EYFP_img(:,:,k),1.5);
    thresh=mean(H2B_filt(:,:,k),'all')+beta*std(H2B_filt(:,:,k),0,'all');
    mask_H2B(:,:,k)=imclearborder(bwareafilt(imbinarize(H2B_filt(:,:,k),thresh),[8e2 5e4]));
    mask_nuc(:,:,k)=imclearborder(bwareafilt(imfill(imbinarize(H2B_filt(:,:,k),thresh),'holes'),[8e2 5e4]));

    area_b(k,1)=sum(mask_nuc(:,:,k),[1 2]).*0.085^2; %um^2

    mat_b(k,1)=sum(EYFP_filt(:,:,k).*mask_nuc(:,:,k),'all');

    local_mask(:,:,k) = bwareafilt(logical(mask_nuc(:,:,k)), [1e3, 5e4]);

    % segment dense phase
    circavgfilt=fspecial('disk',15);
    meanimg=imfilter(EYFP_filt(:,:,k),circavgfilt);
    imgsq=EYFP_filt(:,:,k).^2;
    meansq=imfilter(imgsq,circavgfilt);
    stdevimg=(meansq-meanimg.^2).^0.5;
    threshmat=meanimg+0.5*stdevimg;
    mask_EYFP_final(:,:,k)=bwareafilt(imfill(imbinarize(EYFP_filt(:,:,k),threshmat),'holes'),[80 800]);
    mask_EYFP_final(:,:,k)=imerode(mask_EYFP_final(:,:,k).*local_mask(:,:,k),strel('disk',1));
    mask_EYFP_final(:,:,k)=mask_EYFP_final(:,:,k).*local_mask(:,:,k);
    mask_EYFP_final(:,:,k)=bwareafilt(logical(mask_EYFP_final(:,:,k)),[80 800]);

    % segment dilute phase by dilating the dense phase binary masks
    mask_EYFP_dilate_max(:,:,k) = imdilate(mask_EYFP_final(:,:,k),strel('disk',6));
    mask_EYFP_dilate(:,:,k) = imdilate(mask_EYFP_final(:,:,k),strel('disk',2));
    dilute_mask(:,:,k) = mask_EYFP_dilate_max(:,:,k)-mask_EYFP_dilate(:,:,k);

    % save final mask at each frame as png, can comment this out
    colormap(gray)
    imagesc(EYFP_filt(:,:,k))
    axis equal, axis off
    hold on
    boundary = bwboundaries(mask_EYFP_final(:,:,k));
    for b = 1:length(boundary)
        plot(boundary{b}(:,2), boundary{b}(:,1),'Color','m','LineWidth',1);
        hold on
    end
    dilute = bwboundaries(dilute_mask(:,:,k));
    for c = 1:length(dilute)
        plot(dilute{c}(:,2), dilute{c}(:,1),'Color','b','LineWidth',1);
        hold on
    end
    Image = getframe(gcf);
    imwrite(Image.cdata, sprintf('EYFP_mask_%s_frame%d.png', filename, k));
    close

    % dense and dilute phase intensity
    EYFP_dilute_int(k,1) = mean(nonzeros(dilute_mask(:,:,k).*EYFP_filt(:,:,k)),'all');
    EYFP_dense_int(k,1) = mean(nonzeros(mask_EYFP_dilate(:,:,k).*EYFP_filt(:,:,k)),'all');

    PC(k,1) = EYFP_dense_int(k,1)./EYFP_dilute_int(k,1);

end

end
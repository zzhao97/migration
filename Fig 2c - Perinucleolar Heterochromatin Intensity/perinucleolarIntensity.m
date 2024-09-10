% Confined migration paper JZ 2024
%
%
% - Used in analyzing perinucleolar heterochromatin intensity
% - Reads multi-color tif stacks with H2B and nucleoli channels
% - NPM1 channel is used to segment nucleoli
% - Nucleoli mask are dilated by 6-pixel, then subtracted by nucleoli mask
% to get the perinucleolar heterochromatin region
% - Each frame is saved with perinucleolar heterochromatin boundary
% overlayed
% - Returns average intensity within the perinucleolar heterochromatin mask
%
%
clear; 
close all; 
clc;

output_dir = '/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/2024/20240630_NPM1_perinucleolar_chromatin/perinucleolarIntensityOutputs/';
cd(output_dir)

directory_data='/Users/zhiyuez/Princeton Dropbox/Jessica Zhao/Princeton Brangwynne Lab/Data/2024/20240630_NPM1_perinucleolar_chromatin/whole tif/';

addpath(directory_data)

dir_data = dir([directory_data '*ctrl*.tif']);
fnames_data={};
[fnames_data{1:length(dir_data)}]=dir_data(:).name;

for fnum = 1:length(dir_data) % change dataset
    filename = fnames_data{fnum}; % change dataset
    [pnh_avg_int_temp, area_temp] = analyzePNH(filename);
    Tdata = table(repmat(filename,[length(pnh_avg_int_temp) 1]), pnh_avg_int_temp, area_temp, 'VariableNames',{'Filename','MeanSurfInt','NucArea'});
    writetable(Tdata, sprintf('Perinucleolar_H2B_%s.csv', fnames_data{fnum})) % change filename and dataset
end

% perinucleolar heterochromatin condensation
function [pnh_avg_int, area] = analyzePNH(filename)

% imaging-specific parameters
px = 0.085;
nuc_size_thresh = [7e3 9e4];
NPM1_size_thresh = [400 5500];

% read file and background subtraction
size_img=size(imfinfo(filename),1);
ncolors=3;
nframes=size_img/ncolors;

H2B_img=zeros([size(imread(filename)) nframes]);
NPM1_img=H2B_img;

for k=1:nframes
    H2B_img(:,:,k)=imread(filename,ncolors*k-2);
    NPM1_img(:,:,k)=imread(filename,ncolors*k-1);
end

background_H2B=14;
H2B_img=H2B_img-background_H2B;
H2B_img(H2B_img<0) = 0;

% segment nucleoli
for k=1:nframes

    H2B_img_filt(:,:,k)=imgaussfilt(H2B_img(:,:,k),1);
    thresh_H2B=mean(H2B_img_filt(:,:,k),'all')+0.25*std(H2B_img_filt(:,:,k),0,'all');
    mask_H2B(:,:,k)=imclearborder(bwareafilt(imfill(imbinarize(H2B_img_filt(:,:,k),thresh_H2B),'holes'),nuc_size_thresh));
    area(k,1)=sum(mask_H2B(:,:,k),[1 2]).*px^2; %um^2
    
    NPM1_img_filt(:,:,k)=imgaussfilt(NPM1_img(:,:,k));
    thresh_NPM1=mean(NPM1_img_filt(:,:,k),'all')+2.5*std(NPM1_img_filt(:,:,k),0,'all');
    mask_NPM1(:,:,k)=mask_H2B(:,:,k).*imclearborder(bwareafilt(imfill(imbinarize(NPM1_img_filt(:,:,k),thresh_NPM1),'holes'),NPM1_size_thresh));
    
    % % save mask_NPM1 at each frame as png
    % colormap(jet)
    % imagesc(NPM1_img_filt(:,:,k))
    % axis equal, axis off
    % hold on
    % boundary = bwboundaries(mask_NPM1(:,:,k));
    % for b = 1:length(boundary)
    %     plot(boundary{b}(:,2), boundary{b}(:,1),'Color','m','LineWidth',1);
    %     hold on
    % end
    % Image = getframe(gcf);
    % imwrite(Image.cdata, sprintf('NPM1_seg_%s_frame%d.png', filename, k));
    % close

end

% dilate nucleoli mask by 6-pixel to include perinucleolar heterochromatin
se1 = strel('disk',1);
se2 = strel('disk',6);

for m=1:nframes

    mask_NPM1_dilate(:,:,m) = imdilate(mask_NPM1(:,:,m),se1);
    mask_NPM1_erode(:,:,m) = imerode(mask_NPM1(:,:,m),se2);
    masked_pnh(:,:,m) = mask_NPM1_dilate(:,:,m).*H2B_img_filt(:,:,m);
    masked_surface(:,:,m) = (mask_NPM1_dilate(:,:,m) - mask_NPM1_erode(:,:,m)).*H2B_img_filt(:,:,m);
    
    % % save mask_NPM1 at each frame as png
    % figure
    % colormap(jet)
    % imagesc(H2B_img_filt(:,:,m))
    % axis equal, axis off
    % hold on
    % boundary_dilate = bwboundaries(mask_NPM1_dilate(:,:,m));
    % boundary_erode = bwboundaries(mask_NPM1_erode(:,:,m));
    % for b = 1:length(boundary_dilate)
    %     plot(boundary_dilate{b}(:,2), boundary_dilate{b}(:,1),'Color','m','LineWidth',1);
    %     hold on
    % end
    % 
    % for w = 1:length(boundary_erode)
    %     plot(boundary_erode{w}(:,2), boundary_erode{w}(:,1),'Color','m','LineWidth',1);
    %     hold on
    % end
    % Image = getframe(gcf);
    % imwrite(Image.cdata, sprintf('PNH_seg_%s_frame%d.png', filename, m));  
    % close all

end

% calculate average intensity of surface H2B
masked_surface(masked_surface==0) = NaN;
for j=1:nframes
    pnh_avg_int(j,1) = mean(~isnan(masked_surface(:,:,j)),'all');
end

end
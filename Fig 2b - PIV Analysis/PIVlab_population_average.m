clear;
clc;
close all

directory = '/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/PIVlab data/cell mats/';
% directory = 'E:\Dropbox\Dropbox (Princeton)\Princeton Brangwynne Lab\Data\matlab\PIVlab data\cell mats\';

addpath(directory)
cd(directory)

addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/PIVlab data/')
addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/MatPIV161/postprocessing/')
addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/MatPIV161/src/')
% addpath('E:\Dropbox\Dropbox (Princeton)\Princeton Brangwynne Lab\Data\matlab\MatPIV161\src\')
% addpath('E:\Dropbox\Dropbox (Princeton)\Princeton Brangwynne Lab\Data\matlab\MatPIV161\postprocessing\')

%% Get the population average H2B intensity map
% load images and masks
% img_directory = 'E:\Dropbox\Dropbox (Princeton)\Princeton Brangwynne Lab\Data\matlab\PIVlab data\original sequences\';
% mask_directory = 'E:\Dropbox\Dropbox (Princeton)\Princeton Brangwynne Lab\Data\matlab\PIVlab data\inverted mask sequences\';
img_directory = '/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/PIVlab data/original sequences/';
mask_directory = '/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/PIVlab data/inverted mask sequences/';
addpath(img_directory)
addpath(mask_directory) 

img_dir_mat=dir([img_directory '*.tif']);%find all .mat files in the directory
img_fnames={}; 
[img_fnames{1:length(img_dir_mat)}]=img_dir_mat(:).name; %create cell array containing file names

mask_dir_mat=dir([mask_directory '*.tif']);%find all .mat files in the directory
mask_fnames={}; 
[mask_fnames{1:length(mask_dir_mat)}]=mask_dir_mat(:).name; %create cell array containing file names

cellID = cellfun(@(x) x(1:end-7), img_fnames, 'UniformOutput', false);
unique_arr = unique(cellID);
unique_cellID = cellfun(@(x) strrep(x, '.tif', ''), unique_arr, 'UniformOutput', false);
unique_cellID = unique(unique_cellID);

% calculate average H2B intensity image
H2B_map_mat = zeros([221, 156, length(unique_cellID)]);

for ind = 1:length(unique_cellID)
    current_img_fnames = img_fnames(contains(img_fnames,unique_cellID{ind}));
    current_mask_fnames = mask_fnames(contains(mask_fnames,unique_cellID{ind}));
    img_first_frame = current_img_fnames{1};
    mask_first_frame = current_mask_fnames{1};

    H2B_map_mat(:,:,ind) = getH2BIntFull(img_directory, mask_directory, img_first_frame, mask_first_frame);
end

average_H2B_map_full = mean(H2B_map_mat,3);

% downsize H2B map to match dimension of eta matrix
averageFunction = @(blockStruct) mean(blockStruct.data(:)); % Define a custom function for averaging each block

downsized_H2B_map = blockproc(average_H2B_map_full(4:end-3, 4:end-3), [5,5], averageFunction); % Use blockproc to apply the custom function to each block

%% Load variables from .mat file 
dir_mat=dir([directory '*.mat']);%find all .mat files in the directory
fnames={}; 
[fnames{1:length(dir_mat)}]=dir_mat(:).name; %create cell array containing file names
numfiles=length(dir_mat);

for fnum=1:numfiles 
    filename=fnames{fnum};
    
    x_struct_temp = load(filename,'x');
    x_all_t{fnum} = x_struct_temp.x;
    
    y_struct_temp = load(filename,'y');
    y_all_t{fnum} = y_struct_temp.y;
    
    type_original_temp = load(filename, 'typevector_original');
    type_orig_all_t{fnum} = type_original_temp.typevector_original;
    
    type_filtered_temp = load(filename, 'typevector_filtered');
    type_filt_all_t{fnum} = type_filtered_temp.typevector_filtered;
    
    u_struct_temp = load(filename, 'u_filtered');
    u_all_t{fnum} = u_struct_temp.u_filtered;
    
    v_struct_temp = load(filename, 'v_filtered');
    v_all_t{fnum} = v_struct_temp.v_filtered;
    
    vel_mag_struct_temp = load(filename, 'velocity_magnitude');
    vel_mag_all_t{fnum} = vel_mag_struct_temp.velocity_magnitude;
end

%% Pick out elements passing post-processing criteria based on type struct
for fnum = 1 : numfiles
    filtered_elements = {};
    for tint = 1:size(type_orig_all_t{fnum})
        filtered_elements = type_orig_all_t{fnum}{tint} - (type_filt_all_t{fnum}{tint} - type_orig_all_t{fnum}{tint});
        final_u{fnum}{tint} = u_all_t{fnum}{tint} .* filtered_elements;
        final_v{fnum}{tint} = v_all_t{fnum}{tint} .* filtered_elements;
        final_velmag{fnum}{tint} = vel_mag_all_t{fnum}{tint} .* filtered_elements;
    end
end

%% Calculate time average over 1-min stream for each cell
for fnum = 1 : numfiles
    sumMatrices_u = zeros(size(final_u{fnum}{1}));
    sumMatrices_v = zeros(size(final_v{fnum}{1}));
    sumMatrices_velmag = zeros(size(final_velmag{fnum}{1}));
    
    % Loop through the cell array and accumulate the matrices
    for j = 1:numel(final_u{fnum}) % loop through adjacent frames
        sumMatrices_u = sumMatrices_u + final_u{fnum}{j};
        sumMatrices_v = sumMatrices_v + final_v{fnum}{j};
        sumMatrices_velmag = sumMatrices_velmag + final_velmag{fnum}{j};
    end
    
    % Calculate the mean by dividing by the number of matrices
    mean_u{fnum,1} = sumMatrices_u / numel(final_u{fnum});
    mean_v{fnum,1} = sumMatrices_v / numel(final_v{fnum});
    mean_velmag{fnum,1} = sumMatrices_velmag / numel(final_velmag{fnum});
end

%% Calculate population average of physical parameters
sumAllCell_u = zeros(size(mean_u{1,1}));
sumAllCell_v = zeros(size(mean_v{1,1}));
sumAllCell_velmag = zeros(size(mean_velmag{1,1}));

for fnum = 1 : numfiles
    sumAllCell_u = sumAllCell_u + mean_u{fnum};
    sumAllCell_v = sumAllCell_v + mean_v{fnum};
    sumAllCell_velmag = sumAllCell_velmag + mean_velmag{fnum};
end

all_mean_u = sumAllCell_u / numfiles;
all_mean_v = sumAllCell_v / numfiles;
all_mean_velmag = sumAllCell_velmag / numfiles;

%% coordinate x and y
x_coord = x_all_t{1}{1};
y_coord = y_all_t{1}{1};

%% Calculate shear strain (epsxy) and normal strain (eta) using MatPIV 1.6.1 script
[epsxy,eta]=strain(x_all_t{1}{1},y_all_t{1}{1},all_mean_u,all_mean_v,'centered'); % circulation (3:end-2) and centered (1:end) works well

%% Calculate material derivative of density (product of density and strain)
materialDiffDensity = downsized_H2B_map.*eta;
figure
imagesc(materialDiffDensity), axis equal

%% Visualize physical parameters in color map
addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/cbrewer2/')
addpath('/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/matlab/colorspace/colorspace/')
% addpath('E:\Dropbox\Dropbox (Princeton)\Princeton Brangwynne Lab\Data\matlab\cbrewer2\')
% addpath('E:\Dropbox\Dropbox (Princeton)\Princeton Brangwynne Lab\Data\matlab\colorspace\colorspace')

figure

%h = pcolor(x_coord, y_coord, all_mean_velmag*0.159/0.2); % this is converted to real units
%h = pcolor(x_coord(2:end-1, 2:end-1), y_coord(2:end-1, 2:end-1), epsxy/0.2); % shear strain
h = pcolor(x_coord, y_coord, materialDiffDensity/0.2); % material derivative of density

%colormap(cbrewer2('Greens'))
colormap(cbrewer2('PRGn'))
set(h, 'EdgeColor', 'none');
set(gcf,'Position',[100 100 400 400])
axis off
axis equal

clim([-10 10])
c = colorbar;
c.FontSize = 16;
%c.Label.String = 'Velocity Magnitude (\mum/s)';
c.Color = "black";
c.Location = "north";

%% Visualize flow vector field by Quiver plot (this will be inverted in y)
hold on 
quiver(x_coord, y_coord, all_mean_u, all_mean_v, 2.5, 'color', 'k', 'linewidth', 0.75)
%quiver(x_coord, y_coord, all_mean_u, all_mean_v, 2.5, 'color', 'w', 'linewidth', 0.75)
set(gca, 'color', 'none');
axis off
axis equal

camroll(90)

%% y-axis collaped averaged profile 
% eta_original = eta;
% eta(eta<0) = 0;
% % eta_pos_collapse = sum(eta,2);
% eta_pos_avg = sum(eta,2)./sum(logical(eta),2);
% 
% eta_original(eta_original>0) = 0;
% % eta_neg_collapse = sum(-eta_original,2);
% eta_neg_avg = sum(eta_original,2)./sum(logical(eta_original),2);

mat_diff = materialDiffDensity;
mat_diff_copy = materialDiffDensity;
mat_diff(mat_diff<0) = 0;
diff_pos_avg = sum(mat_diff,2)./sum(logical(mat_diff),2);

mat_diff_copy(mat_diff_copy>0) = 0;
diff_neg_avg = sum(mat_diff_copy,2)./sum(logical(mat_diff_copy),2);

figure
plot(linspace(0,length(mat_diff)*0.159*5,length(mat_diff)),flip(diff_pos_avg),'color','#5aae61','linewidth',1.5)
hold on
plot(linspace(0,length(mat_diff)*0.159*5,length(mat_diff)),flip(diff_neg_avg),'color','#9970ab','linewidth',1.5)
set(gca,'FontSize',25)
set(gca,'LineWidth',1)
xlim([0 33])
xlabel('position along x-direction (\mum)')
ylabel('Material Derivative of Density')
set(gcf,'Position',[400 400 600 400])
legend('positive','negative')

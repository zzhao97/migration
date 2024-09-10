%batch process
clear
clc
close

directory='/Users/zhiyuez/Dropbox (Princeton)/Princeton Brangwynne Lab/Data/20231211_DZNep_24hours_corelets/488Conversion/';
addpath(directory)
cd(directory)

dirobj=dir([directory '*.tif']);%find all tifs in the directory
fnames={};
[fnames{1:length(dirobj)}]=dirobj(:).name; %create cell array containing file names
numfiles=length(dirobj);

for fnum=1:numfiles %loops through file names
    filename=fnames{fnum};
    [Ch1Int_temp, Ch2Int_temp]= get488Intensity(directory, filename); %apply function to each image
    Ch1Int{fnum,1}=Ch1Int_temp; %write outputs
    Ch2Int{fnum,1}=Ch2Int_temp;
end
%%
%savedata
Tdata=table(cell2mat(Ch1Int),cell2mat(Ch2Int),'VariableNames',{'Channel1','Channel2'}) %generate output table
writetable(Tdata,[directory 'output_data.csv'])
%% fitting and plot
degree = 1;
coefficients = polyfit(cell2mat(Ch1Int), cell2mat(Ch2Int), degree);
x = cell2mat(Ch1Int);
y = cell2mat(Ch2Int);
xfit = linspace(min(x), max(x), 100);

figure %plot the results
plot(x,y,'or', xfit, polyval(coefficients, xfit), '-b')
xlabel('I_{1}')
ylabel('I_{2}')
set(gca,'FontSize',25)
axis square
ylim([0 250])
%%
function [Ch1Int, Ch2Int]= get488Intensity(directory, filename)
fullfilename=[directory filename];
backlvl=100;
ch2=imread(fullfilename,2);
ch1=imread(fullfilename,1);
ch2 = double(ch2)-backlvl;
ch1 = double(ch1)-backlvl;

ch2(ch2<0) = 0;
ch1(ch1<0) = 0;

ch1_filt = imgaussfilt(ch1);
ch2_filt = imgaussfilt(ch2);

meanimg=imfilter(ch1,fspecial('disk',60));
imgsq=ch1_filt.^2;
meansq=imfilter(imgsq,fspecial('disk',60));
stdevimg=(meansq-meanimg.^2).^0.5;
threshmat=meanimg+0.25.*stdevimg;
ch1_mask_core=bwareafilt(imbinarize(ch1_filt,threshmat),[9e3 9e4]);

Ch1Int = mean(ch1.*ch1_mask_core,'all');
Ch2Int = mean(ch2.*ch1_mask_core,'all');
end
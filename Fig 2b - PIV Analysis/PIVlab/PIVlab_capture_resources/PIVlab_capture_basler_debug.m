
% Prepare camera
delete(imaqfind); %clears all previous videoinputs
hwinf = imaqhwinfo;
info = imaqhwinfo(hwinf.InstalledAdaptors{1});
basler_name = info.DeviceInfo.DeviceName;
basler_supported_formats = info.DeviceInfo.SupportedFormats;
basler_vid = videoinput(info.AdaptorName);
basler_settings = get(basler_vid);
basler_settings.Source.DeviceLinkThroughputLimitMode = 'off';

%basler_vid.ROIPosition=[700,300,640,480]


%% set camera parameters for free run or triggered acquisition

triggerconfig(basler_vid, 'hardware');
basler_settings.TriggerSource = 'Line3';
basler_settings.Source.ExposureMode = 'TriggerWidth';
basler_settings.Source.TriggerSource ='Line3';
basler_settings.Source.TriggerSelector='FrameStart';
basler_settings.Source.TriggerMode ='On';
basler_settings.Source.ExposureOverlapTimeMax = basler_settings.Source.SensorReadoutTime;

image_handle_basler=imagesc(zeros(basler_settings.VideoResolution(2),basler_settings.VideoResolution(1)),[0 2^8]);
axis image
axis off


basler_frames_to_capture = 100;
basler_vid.FramesPerTrigger = basler_frames_to_capture;
flushdata(basler_vid);
start(basler_vid);
preview(basler_vid,image_handle_basler)
drawnow;

serpo=serialport("COM1",9600);
configureTerminator(serpo,"CR/LF");
writeline(serpo,"TALKINGTO:SYNC_UOA1;FREQ:168;CAM:0;ENER:0;ener%:100;F1EXP:0;INTERF:250;EXTDLY:0;EXTSKP:0;LASER:enable{13}");

while basler_vid.FramesAcquired < (basler_frames_to_capture)
	drawnow limitrate
end
stoppreview(basler_vid)
stop(basler_vid);


basler_data = getdata(basler_vid); %ruft alle Frames in RAM ab. Frame 1,2,3 sind müll

writeline(serpo,"TALKINGTO:SYNC_UOA1;FREQ:168;CAM:0;ENER:0;ener%:100;F1EXP:0;INTERF:250;EXTDLY:0;EXTSKP:0;LASER:disable{13}");
serpo=[];
clear serpo


ImagePath ='C:\Users\trash\Desktop\PIV_DATA\basler';
cntr=0;
for image_save_number=1:2:size(basler_data,4)-1
	imgA_path=fullfile(ImagePath,['PIVlab_' sprintf('%4.4d',cntr) '_A.tif']);
	imgB_path=fullfile(ImagePath,['PIVlab_' sprintf('%4.4d',cntr) '_B.tif']);
	imwrite(basler_data(:,:,:,image_save_number),imgA_path,'compression','none'); %tif file saving seems to be the fastest method for saving data...
	imwrite(basler_data(:,:,:,image_save_number+1),imgB_path,'compression','none');
	cntr=cntr+1;
end
close all
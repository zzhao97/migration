
% Prepare camera
delete(imaqfind); %clears all previous videoinputs
clear
hwinf = imaqhwinfo;
info = imaqhwinfo(hwinf.InstalledAdaptors{1});
flir_name = info.DeviceInfo.DeviceName;
flir_supported_formats = info.DeviceInfo.SupportedFormats;
flir_vid = videoinput(info.AdaptorName);
flir_settings = get(flir_vid);
flir_settings.Source.DeviceLinkThroughputLimit = flir_settings.Source.DeviceLinkSpeed;


%% set camera parameters for free run or triggered acquisition

flir_settings.Source.TriggerSource='Line2';
flir_settings.Source.TriggerActivation='RisingEdge';
flir_settings.Source.TriggerMode='On';
flir_settings.Source.TriggerSelector='FrameStart';
flir_settings.Source.ExposureMode = 'Timed';
triggerconfig(flir_vid, 'hardware');
%trigger width setting does not work on the FLIR, therefore using exposure time.
%min blind time between frames is 400 µs
debug_frame_rate=60;
flir_settings.Source.ExposureTime = floor(1/debug_frame_rate*1000*1000-405);



%%set line3 to output exposureactive signal
flir_settings.Source.LineSelector='Line3';
flir_settings.Source.LineSource = 'ExposureActive';
flir_settings.Source.LineMode = 'Output';
flir_settings.Source.LineInverter='False';



image_handle_flir=imagesc(zeros(flir_settings.VideoResolution(2),flir_settings.VideoResolution(1)),[0 2^8]);
axis image
axis off


flir_frames_to_capture = 100;
flir_vid.FramesPerTrigger = flir_frames_to_capture;
flushdata(flir_vid);
start(flir_vid);
preview(flir_vid,image_handle_flir)
drawnow;

%{
serpo=serialport("COM1",9600);
configureTerminator(serpo,"CR/LF");
writeline(serpo,"TALKINGTO:SYNC_UOA1;FREQ:168;CAM:0;ENER:0;ener%:100;F1EXP:0;INTERF:250;EXTDLY:0;EXTSKP:0;LASER:enable{13}");
%}
while flir_vid.FramesAcquired < (flir_frames_to_capture)
	drawnow limitrate
end
stoppreview(flir_vid)
stop(flir_vid);


flir_data = getdata(flir_vid); %ruft alle Frames in RAM ab. Frame 1,2,3 sind müll

%{
writeline(serpo,"TALKINGTO:SYNC_UOA1;FREQ:168;CAM:0;ENER:0;ener%:100;F1EXP:0;INTERF:250;EXTDLY:0;EXTSKP:0;LASER:disable{13}");
serpo=[];
clear serpo
%}

ImagePath ='C:\Users\trash\Desktop\PIV_DATA\flir';
cntr=0;
for image_save_number=1:2:size(flir_data,4)-1
	imgA_path=fullfile(ImagePath,['PIVlab_' sprintf('%4.4d',cntr) '_A.tif']);
	imgB_path=fullfile(ImagePath,['PIVlab_' sprintf('%4.4d',cntr) '_B.tif']);
	imwrite(flir_data(:,:,:,image_save_number),imgA_path,'compression','none'); %tif file saving seems to be the fastest method for saving data...
	imwrite(flir_data(:,:,:,image_save_number+1),imgB_path,'compression','none');
	cntr=cntr+1;
end
close all
input = getDirectory("input directory");
original_output = getDirectory("output directory");
mask_output = getDirectory("mask directory")

action(input);

function action(input){
	fileList = getFileList(input)
	for (i = 0; i < fileList.length; i++){
		open(input + fileList[i]);
		//run("Image Sequence... ", "format=TIFF name=" + fileList[i] + " start=1 digits=3 save=["+original_output+"]");
		selectWindow(fileList[i]);
		run("Duplicate...", "duplicate");
		run("Smooth", "stack");
		setAutoThreshold("Default dark");
		//run("Threshold...");
		setAutoThreshold("IsoData dark");
		setOption("BlackBackground", true);
		run("Convert to Mask", "method=IsoData background=Dark calculate black");
		run("Fill Holes", "stack");
		run("Erode", "stack"); run("Erode", "stack");
		run("Analyze Particles...", "size=75-Infinity show=Masks display exclude stack");
		run("Invert", "stack");
		run("Image Sequence... ", "format=TIFF name=InvertMask" + fileList[i] + " start=1 digits=3 save=["+mask_output+"]");
		while(nImages>0)close();
	}
}

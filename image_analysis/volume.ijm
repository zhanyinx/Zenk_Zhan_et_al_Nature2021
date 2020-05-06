
macro "volumns 1 and 2 [F10]"{

	home=getDirectory("home")
	name1=getTitle();

	//Set standard ROI xy
	rectangle_size=60
	getSelectionBounds(x, y, width, height);
	diff_width = (rectangle_size-width)/2
	diff_height = (rectangle_size-height)/2
	makeRectangle(x-diff_width, y-diff_height, rectangle_size, rectangle_size);

	//Set standard zstacks
	Stack.getPosition(channel, slice, frame);
	getDimensions(w, h, channels, slices, frames);
	low = maxOf(0,slice-3);
	high = minOf(slices,slice+3);
  	run("Duplicate...", "duplicate channels=1-2 slices="+low+"-"+high);
	
	//Make binary
	run("Make Binary", "method=MaxEntropy background=Dark calculate black");

	
	//Split channels
  	name=getTitle;
  	run("Split Channels");
  	img1 = "C1-"+name;
  	img2 = "C2-"+name;

	
	//Distance between cm
	selectWindow(img1);
	vc1 = 0.
	vc2 = 0.
	getDimensions(w, h, channels, slices, frames);
	
	for (i = 0; i < slices; i++) {
		run("Clear Results");
		
		selectWindow(img1);
		Stack.setSlice(i)
  		run("Set Measurements...", "integrated nan redirect=None decimal=3");
		run("Measure");
		
		selectWindow(img2);
		Stack.setSlice(i)
  		run("Set Measurements...", "integrated nan redirect=None decimal=3");
		run("Measure");
		
		
		tmpvc1 = getResult("RawIntDen",0);
		tmpvc2 = getResult("RawIntDen",1);

		vc1 = vc1 + tmpvc1 / 255;
		vc2 = vc2 + tmpvc2 / 255;

	}
	File.append(vc1+ " " + vc2 + " " + (x-diff_width)+"_" + y-diff_height+ "_" + low+ "_" + high ,home+"/"+name1+".volumes.txt")
	
	//Cleaning
	selectWindow(img1);
	close();
	selectWindow(img2);
	close();
	
}
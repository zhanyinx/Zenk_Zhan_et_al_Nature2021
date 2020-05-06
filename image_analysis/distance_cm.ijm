macro "distance cm channels 1 and 2 [F5]"{

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
	
	//Transform to binary
	run("Make Binary", "method=Default background=Dark calculate black");

	
	//Split channels
  	name=getTitle;
  	run("Split Channels");
  	img1 = "C1-"+name;
  	img2 = "C2-"+name;

	
	//Distance between cm
	selectWindow(img1);
	dx = 0.
	dy = 0.
	counter = 0
	getDimensions(w, h, channels, slices, frames);
	
	for (i = 0; i < slices; i++) {
		run("Clear Results");
		
		selectWindow(img1);
		Stack.setSlice(i)
  		run("Set Measurements...", "area standard centroid center area_fraction stack nan redirect=None decimal=3");
		run("Measure");
		
		selectWindow(img2);
		Stack.setSlice(i)
  		run("Set Measurements...", "area standard centroid center area_fraction stack nan redirect=None decimal=3");
		run("Measure");
		
		
		tmpxm1 = getResult("XM",0);
		tmpym1 = getResult("YM",0);
		tmpxm2 = getResult("XM",1);
		tmpym2 = getResult("YM",1);

		dx = dx + pow(tmpxm1-tmpxm2,2);
		dy = dy + pow(tmpym1-tmpym2,2);
		counter = counter + 1;

	}
	d = sqrt((dx + dy)/counter)
	File.append(d+ " "+ (x-diff_width)+"_" + y-diff_height+ "_" + low+ "_" + high ,home+"/"+name1+".txt")
	
	//Cleaning
	selectWindow(img1);
	close();
	selectWindow(img2);
	close();
	
}
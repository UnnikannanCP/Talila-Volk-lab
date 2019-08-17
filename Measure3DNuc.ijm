/* 
 *  Measure3DNuc - Segment 3D nucleus ferom a selected chanel and Quantify shape and Intensity of other channels  
 *  For: Daria Amiad-Pavlov, Talila Volk
 *  By: Ofra Golani, january 2019
 *  
 *  v0.1, July 2019: fix problem with multiple run (3DManager multiple instances) 
 *  
 *  Usage Instructions
 * ====================
 * 1 Set Parameters as needed
 * 
 * 2 Click the Run Button, The macro will prompt you to select an input stack for analysis 
 * 
 * 3 The macro will prompt you to check Threshold validity, 
 *   the easiest way to do it is to put the "Nuc" image next to the "composite" and "MaxProject" images, 
 *   synchronize the windows using Anlyze=>Tools=>Synchronize Windows and scroll through the slices to see the slices. 
 *   Click OK when done
 *   If you are not satisfied try to select another Threshold (Image=>Adjust=>Threshold... , global using all the stack, stop the macro, and set NucAutoThMethod to the new method
 *   Note however, that it is recommended to use the same method (or same fixed threshold,not implemented here though) for all the files in the experiment  
 *   
 * 4 The macro will add the segmented Nuclei to the 3D ROI Manager and stop to let you browse the resulted objects and discard unneeded ones.  
 *   the labeled Image "NucProjectLabelTmp" can help you with this
 *   Nuc that are not fully included in the stack should be discraded, as well as those that are not in the correct muscle.  
 *   To discard a Nuc, you need to Click LiveROI in the 3D Roi Manager, go through individual ROIs and ERASE the unneded one. 
 *   You can display on other images, and scroll through the images, BUT MAKE SURE to ERASE when NucObjectsMap image is active !
 *   Note that after each Erase you'll need to activate the LiveROI again  
 *   You can inspect Individual objects in 3D by clicking on "3D Viewer" from the 3D Manager. You'll need to use Edit=>Change color to see it better. This does not work for the whole image
 *   When Done click OK 
 *   
 * 5 The macro will stop again to let you choose tables and features to color code some the Nuc based on their intensity content.
 *   this time the instructions will appear, and you'll be prompted for the actual action after you click OK, 
 *   Set Min/Max as desired, but it is advised to select smaller Min value than automatically suggested 
 *   
 * 6 At the end you'll "Done !" is printed to the Log window
 * 
 * Note: You can ignore the console with Red messages if it appears 
 * 
 *  General settings
 * -------------------------------
 * Channel Settings: the macro is built for to deal with flexible channel setting. 
 * It is assumed that there is at least one channel for segmentation and optionally channel for DNA normalization. 
 * All other channels are used for Intensity quantification 
 * 
 * NucSegChannel - the number of the channel to use for segmentation, usually this should be the Lamin channel. channel numbers start from 1
 * NucDNAChannel - the numbe of the channel to be used for DNA content normalization. usually this should be the DAPI channel. channel numbers start from 1, if set to -1, no normalization is done
 * 
 * Segmentation Parameters 
 * -------------------------------
 * RollingBallRadius 	- radius in um of the 2D Roling ball used for background subtraction (BGS) of all channels. note that it should be larger than the radius of the Nuclei
 * NucBlurSig			- sigma used for Gaussian Blur of the Nuc signal after BGS and before threshold. in Pixels 
 * NucAutoThMethod 		- name of the method used for setting global threshold for the whole stack of the processed Nuc , default is Utsu
 * MinNucSize 			- minimal size in um^3 of Nuc, all smaller objects are discarded 
 * MaxNucSize 			- maximal size in um^3 of Nuc, all larger objects are discarded 
 * ExcludeNucOnEdge 	- either 0 or 1, 1 indicate to discard all the objects that touchany of the stack borders: in x y or z 
 * 
 * Quantification Parameters
 * -------------------------------
 * RollingBallRadius - the same parameter mentioned above, influence both segmentation and quantification. Background subtraction is done using the same RollinBallradius for all the channels
 * 					   and total intensity is quantified before and after background subtraction
 * 					   
 * Operational Parameters
 * -------------------------------
 * fileExtension 					- input file extension 
 * ResultsSubFolder 				- name of sub folder of the original location, in which the results will be saved (default is "Results")
 * CheckThresholdValidity 			- either 0 or 1 (default), 1 means wait and let the user inspect the thresholded objects compared to the data, to verify segmentation is correct;
 * ManualEditNucsegmentation 		- either 0 or 1 (default), 1 means stop to let the user erase some of the Nuc
 * Save3DRoiMangerSeparateTables 	- either 0 or 1 
 *
 * Output Files
 * -------------------------------
 * All output files are saved in a subfolder under the original location.  
 * All result files are start with the original name and have suffix based on the content.
 * Suppose the file name is XX
 * 
 * - XX_Measure3DNucPrms.txt - 
 * - XX_NucObjectsMap.tif  				- stack of labeled Nuc, each Nuc gets different number, color code is arbitrary (galesby inverted LUT)
 * - XX_NucObjectsMapProjectLabel.tif 	- a MaxProject image with object labels, this is usefull to identify objects
 * - XX_NucRoi3D.zip  					- the 3D ROIs that correspond to the segmented Nuc This File can be opened in the 3D Roi Manger for inspection and verification of quantification
 * - XX_NucSummaryTable.csv - result table with one line for each Nuc, include the most important measurements: 
 *                            Centroid (pixels), Volume (um^3), number of Voxels, Surface Area, elongation 1 & 2, 
 *                            Total intensity of each quantified channel after BGS optonally normalized by teh DNA channel
 * 
 * If Save3DRoiMangerSeparateTables=1
 * - M_XX_Measure3D.csv 		- table with all shape measurements from 3D Manager, one line per object, Tab sepearted (In excel,  select all, Data=>Txt-to-columns choose Tab delimeter)
 * - Q_XX_Cn_Bgs_Quantif3D.csv 	- tables for each channel with all intensity measurements from 3D manager, one line per object, Tab sepearted 
 * 	
 * 	If SaveMorphoLibJSeparateTables=1 (However it is redundant when you save the 3D Manager tables)
 * 	XX_Cn_Bgs-intensity-measurements.csv tables for each channel with all intensity measurements from MorphoLibJ, one line per object
 * 	
 * 	XX_ColorCodeObjects_Y...		- color coded (z-project) objects, using the Y feature as color code 
 * 	
 * Results Inspection
 * ====================
 * - Open the original file 
 * - Plugins=>3D=>3D Manager, use the Open button to open XX_NucRoi3D.zip, Click Live ROI
 * - You can get all the voxels (and values) of a given object by selecting the desired channel and clicking "List Voxels" (their sum is the the total Intensity Orig of this channel in the summary table) 
 * 
 *  Workflow
 * ====================
 * - Segment the 3D Nuclei 
 * 		- select proper Channel based on setting of NucSegChannel
 * 		- background subtraction using 2D RollingBall (RollingBallRadius)
 * 		- Gaussian Blur (NucBlurSig)
 * 		- Apply Auto global Threshold for the whole stack (NucAutoThMethod) - At this point let the user 
 * 		- Connected componnet analysis , discard objects based on size [MinNucSize-MaxNucSize]
 * 		- Add Object to 3D Roi Manager
 * 		- Let the user Erase some of the objects 
 * 		- Save the objects
 * - Quantify Shape and intensity of each Nuclei 
 * 		- Quantify 3D shape of Nuc
 * 		- Quantify total Intensity of the DNA Channel before abd after 2D background subtraction (optional)
 * 	 	- Quantify total intensity of each of the Quantify channels before and after 2D background subtraction 
 * 	 	- Optionally normalize (divide) the total intensity of each channel after BGS by the total intensity of the DNA channel 
 * 	 	- Record all the measurements in a summary table with one line for each Nuc
 * 	 	- for convinience additional measurements are saved in separate tables 
 *  
 *  
 *  Dependencies: should be cited when publishing using the macro
 *  0) Fiji. 
 *     Ref: https://imagej.net/Citing 
 *     Citation: Schindelin, J.; Arganda-Carreras, I. & Frise, E. et al. (2012), "Fiji: an open-source platform for biological-image analysis", Nature methods 9(7): 676-682, PMID 22743772, doi:10.1038/nmeth.2019 (on Google Scholar). 
 *  1) 3D Object Counter plugin (Installed by default as part of Fiji)    
 *     Ref: https://imagej.net/3D_Objects_Counter
 *     Citation: S. Bolte & F. P. Cordelières, A guided tour into subcellular colocalization analysis in light microscopy, Journal of Microscopy, Volume 224, Issue 3: 213-232 
 *  2) 3D ImageJ Suite (Installed by default with Fiji, 3D ImageJ Suite should be selected in Update sites). 
 *     Ref: http://imagejdocu.tudor.lu/doku.php?id=plugin:stacks:3d_ij_suite:start 
 *     Citation: J. Ollion, J. Cochennec, F. Loll, C. Escudé, T. Boudier. (2013) TANGO: A Generic Tool for High-throughput 3D Image Analysis for Studying Nuclear Organization. Bioinformatics 2013 Jul 15;29(14):1840-1. doi: http://dx.doi.org/10.1093/bioinformatics/btt276  
 *  3) MorphoLibJ (select IJPB-Plugins from update sites), 
 *     Ref: https://imagej.net/MorphoLibJ
 *     Citation: Legland, D.; Arganda-Carreras, I. & Andrey, P. (2016), "MorphoLibJ: integrated library and plugins for mathematical morphology with ImageJ", Bioinformatics (Oxford Univ Press) 32(22): 3532-3534, PMID 27412086, doi:10.1093/bioinformatics/btw413 (on Google Scholar). 
 */

// =================== Parameters =====================================================
// Quantification Parameters
//-------------------------
// All other channels are used for Quantification
var NucSegChannel = 1; // Channel for segmentation, usualy Lamin channel
var NucDNAChannel = 4; // Channel for total intensity normalization, -1 means don't normalize
var ChColor = newArray("Blue", "Red", "Green", "Grays");

// Nuc segmentation parameters
var RollingBallRadius = 20; //um
var NucBlurSig=2; 			// pixels
var NucAutoThMethod = "Otsu";
var MinNucSize = 100; 		// um^3
var MaxNucSize = 10000;		// um^3
var ExcludeNucOnEdge = 0;

// Operational Parameters
//-------------------------
var fileExtension = ".czi";
var ResultsSubFolder = "Results";
var NucSummaryTableName = "NucSummaryTable";
var CheckThresholdValidity = 1; 
var ManualEditNucsegmentation = 1;
var SaveMaxProject = 1;
var Save3DRoiMangerSeparateTables = 1;
var ColorCodeNormIntensity = 1;

var MorphoLibJIntMeasure_DNAChannel = 0;	
var MorphoLibJIntMeasure_QuantChannel = 0;
var SaveMorphoLibJSeparateTables = 0; 		// No need, redundant to 3D manager Tables
var Label3DFlag = 1;
var OpenIn3DViewer = 0; //this does not work currently, keep it 0
var DebugFlag = 0;
var CleanupFlag = 1;

//==================== Main Code don't change beyond this point =======================

// Helper Global variables
//=====================================================================================
var QuantChannels;
var nQuantChannels;

// Initialization
//=====================================================================================
print("\\Clear"); // clear the log window
if (isOpen(NucSummaryTableName))
{
	selectWindow(NucSummaryTableName);
	run("Close");
}
//call("ij3d.ImageJ3DViewer.close");

// Select the image, create output subfolder
//=====================================================================================
file_name=File.openDialog("Please select an image to analyze");
//open(file_name);
//print(file_name);
//print(File.getParent(file_name)); 
run("Bio-Formats Importer", "open=["+file_name+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT stitch_tiles");
directory = File.getParent(file_name); 
print(directory);
resFolder = directory + File.separator + ResultsSubFolder + File.separator; 
File.makeDirectory(resFolder);

// Refer to original image
orig=getTitle();
origNameNoExt = replace(orig, fileExtension, "");

// save Parameters
SavePrms(resFolder, origNameNoExt);


// Setup Quantification Channels and Prepare independent Channel Images & channels after background subtraction (bgs) for quantification
// Channel Images and bgs images are called: Cn & Cn_bgs respectively
//=====================================================================================
SetupChannelImages(orig); 

// Nuc Segmentation: blur, threshold (2D, stack histogram), 3D connected componnet analysis, add to 3D Roi Manager , output is "NucObjectsMap"
//=====================================================================================
NucSegmentation();

// Manual Editing of Nuc segmentation
//=====================================================================================
selectWindow("NucObjectsMap");
setSlice(round(nSlices/2));
if (ManualEditNucsegmentation == 1)
{
	if (Label3DFlag == 1)
	{
		selectWindow("NucObjectsMap");
		run("Select None");
		run("Duplicate...", "title=NucLabelTmp duplicate");
		Ext.Manager3D_Label();
		run("Z Project...", "projection=[Max Intensity]");
		rename("NucProjectLabelTmp");
	}
	selectWindow("NucObjectsMap");
	waitForUser("Make sure Live Roi is always On, Click on Label\n Erase non-valid Nuceli, using 'Erase' button\nMake sure the NucObjectsMap is selected, you can compare to the MaxProject image\nWhen done click 'OK'");
}
run("Select None");
Ext.Manager3D_Count(nNuc);
print("Number of Nuc after manual editing is: ",nNuc);

// Save 3D Objects after segmentation
//selectWindow("composite");
Ext.Manager3D_SelectAll();
Ext.Manager3D_Save(resFolder+origNameNoExt+"_NucRoi3D.zip");

// Measurments
// Intensity and Shape Quantification from 3D Manager
// Intensity is normalized by DNA channel 
//=====================================================================================

// if list is not visible please refresh list by using Deselect
Ext.Manager3D_Select(0);
Ext.Manager3D_DeselectAll();

// number of results, and arrays to store results
//Ext.Manager3D_Count(nNuc);
// get object labels
label=newArray(nNuc);
cx=newArray(nNuc);
cy=newArray(nNuc);
cz=newArray(nNuc);
vol=newArray(nNuc);
surf=newArray(nNuc);
nvox=newArray(nNuc);
elon1=newArray(nNuc);
elon2=newArray(nNuc);
totDNAInt=newArray(nNuc);

selectWindow("NucObjectsMap");
// loop over objects - location & shape quantification
for(i=0;i<nNuc;i++){
	Ext.Manager3D_GetName(i, label[i]);
	Ext.Manager3D_Centroid3D(i,cx[i],cy[i],cz[i]);	// Centroid
	Ext.Manager3D_Measure3D(i,"Vol",vol[i]); 		// volume
	Ext.Manager3D_Measure3D(i,"Surf",surf[i]); 		// surface
	Ext.Manager3D_Measure3D(i,"NbVox",nvox[i]); 	// number of voxels
	Ext.Manager3D_Measure3D(i,"Elon1",elon1[i]); 	// elongation 1
	Ext.Manager3D_Measure3D(i,"Elon2",elon2[i]); 	// elongation 2
}
// Show measured results in "Results" table
Array.show("Results",label, cx, cy, cz, vol, surf, nvox, elon1, elon2);

// Intensity Quantification: loop over objects to print individual values into tables 
// to support different number of channels and optional normalization
for(i=0;i<nNuc;i++)
{
	if (NucDNAChannel != -1) 
	{
		ChName = "C"+NucDNAChannel;
		selectWindow(ChName);
		run("Select None");
		Ext.Manager3D_Quantif3D(i,"IntDen",quantif);
		//totDNAInt[i] = quantif;
		setResult("totDNAIntOrig", i, quantif);
		
		bgsName = "C"+NucDNAChannel+"_Bgs";
		selectWindow(bgsName);
		run("Select None");
		Ext.Manager3D_Quantif3D(i,"IntDen",quantif);
		totDNAInt[i] = quantif;
		setResult("totDNAIntBgs", i, quantif);
	}
	for (n = 0; n < nQuantChannels; n++)
	{
		id = QuantChannels[n];
		// measure and save total original intensity
		ChName = "C"+id;
		selectWindow(ChName);
		run("Select None");
		Ext.Manager3D_Quantif3D(i,"IntDen",quantif);
		ColName = "totInt"+ChName+"Orig";
		setResult(ColName, i, quantif);

		// Measure and save total intensity after background subtraction
		bgsName = "C"+id+"_Bgs";
		selectWindow(bgsName);
		run("Select None");
		Ext.Manager3D_Quantif3D(i,"IntDen",quantif);
		ColName = "totInt"+bgsName;
		setResult(ColName, i, quantif);
		if (NucDNAChannel != -1) 
		{
			ColName = "totInt"+bgsName+"Norm";
			NormVal = quantif / totDNAInt[i];
			setResult(ColName, i, NormVal);
		}
	}
}

// rename table to NucSummaryTableName
Table.rename("Results", NucSummaryTableName);
saveTable(resFolder, origNameNoExt, NucSummaryTableName, 0);

// For comparison - save individual quantification tables, later on this part of the code can be deleted 
// Measure Intensity, start with DNA normalization channel if exist, then go to Quantif Channels
//=====================================================================================
if (Save3DRoiMangerSeparateTables)
{
	selectWindow("C"+NucSegChannel+"_Bgs");
	// Measure shape 
	Ext.Manager3D_Measure();
	Ext.Manager3D_SaveResult("M",resFolder+origNameNoExt+"_Measure3D.csv");
	if (DebugFlag == 0) Ext.Manager3D_CloseResult("M");
	//waitForUser("After saving Measure3D");
	
	// Measure Intensity
	if (NucDNAChannel > 0) 
	{
		ChName = "C"+NucDNAChannel+"_Bgs";
		MeasureAndSave3DIntensity(resFolder, origNameNoExt, ChName);
	}
	for (n = 0; n < nQuantChannels; n++)
	{
		ChName = "C"+QuantChannels[n]+"_Bgs";
		MeasureAndSave3DIntensity(resFolder, origNameNoExt, ChName);
	}
	//waitForUser("After saving Quantif3D");
}

// Color code nuc by norm total intensity of DNA channel and save table
if ((ColorCodeNormIntensity == 1) && (NucDNAChannel > -1))
{
	ColorCodeFromTable("totDNAIntBgs", NucSummaryTableName, 0, 0);
	for (n = 0; n < nQuantChannels; n++)
	{
		ColorCodeFromTable("totInt"+"C"+QuantChannels[n]+"_Bgs"+"Norm", NucSummaryTableName, 0, 0);
	}
}


// Intensity and Shape Quantification from MorphoLibJ 
// may be redundant, usefull only if color by Mean Intensity is needed
//=====================================================================================
if (SaveMorphoLibJSeparateTables)
{
	selectWindow("NucObjectsMap");
	run("Select None");
	run("Region Morphometry");
	saveTable(resFolder, origNameNoExt, "NucObjectsMap-Morphometry", 1);
	//waitForUser("After save MorphoLibJ Morphometry");
	
	// Color code nuc by intensity mean of DNA channel and save table
	if ((MorphoLibJIntMeasure_DNAChannel == 1) && (NucDNAChannel > -1))
	{
		MorphoLibJIntMeasureAndColorCode(resFolder, origNameNoExt, NucDNAChannel, "Mean");
	}
	
	if (MorphoLibJIntMeasure_QuantChannel == 1)
	{
		for (n = 0; n < nQuantChannels; n++)
		{
			MorphoLibJIntMeasureAndColorCode(resFolder, origNameNoExt, QuantChannels[n], "Mean");
		}
	}
	//waitForUser("After save MorphoLibJ Intensity Measure");
}

// Save Labeled images for Quality Control 
//=====================================================================================
// Projection + label numbers of final Nuc
selectWindow("NucObjectsMap");
run("Select None");
run("Duplicate...", "title=NucObjectsMap_ForLabel duplicate");
Ext.Manager3D_Label();
run("Z Project...", "projection=[Max Intensity]");
saveAs("Tiff", resFolder+origNameNoExt+"_NucObjectsMapProjectLabel.tif");

// save Nuc Object map stack
selectWindow("NucObjectsMap");
saveAs("Tiff", resFolder+origNameNoExt+"_NucObjectsMap.tif");
rename("NucObjectsMap");

// Display in 3D Viewer - this does not work currently, keep OpenIn3DViewer=0
if (OpenIn3DViewer)
{
	selectWindow("NucObjectsMap");
	run("3D Viewer");
	call("ij3d.ImageJ3DViewer.setCoordinateSystem", "false");
	call("ij3d.ImageJ3DViewer.add", "NucObjectsMap", "None", "NucObjectsMap", "0", "true", "true", "true", "1", "0");
	call("ij3d.ImageJ3DViewer.add", "composite",     "None", "composite",     "0", "true", "true", "true", "1", "0");
	call("ij3d.ImageJ3DViewer.resetView");
		
	call("ij3d.ImageJ3DViewer.record360");
	selectWindow("Movie");
	run("AVI... ", "compression=JPEG frame=10 save="+resFolder+origNameNoExt+"_NucObjects3DMovie.avi");
	close();	
}
// Cleanup: close all images, close all tables 
if (CleanupFlag == 1)
{
	run("Close All");
	Ext.Manager3D_CloseResult("M");
	Ext.Manager3D_CloseResult("Q");
	Ext.Manager3D_Reset();
	Ext.Manager3D_Close();
	selectWindow(NucSummaryTableName);
	run("Close");
}
print("Done !");

//==================== End of Main Code =====================================================
//==================== Helper Functions =====================================================

// Measure 3D Intensity from 3D Suite, save and close the table
//---------------------------------------------------------
function MeasureAndSave3DIntensity(resFolder, origNameNoExt, ChName)
{
	selectWindow(ChName);
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveResult("Q",resFolder+origNameNoExt+"_"+ChName+"_Quantif3D.csv");
	if (DebugFlag == 0)
		Ext.Manager3D_CloseResult("Q");
}
	
// Measure Intensity using MorphoLibJ plugin, Color code nuc by intensity mean of N channel and save table
//---------------------------------------------------------
function MorphoLibJIntMeasureAndColorCode(resFolder, origNameNoExt, Channel, Feature)
{
	bgsName = "C"+Channel+"_Bgs";
	run("Intensity Measurements 2D/3D", "input="+bgsName+" labels=NucObjectsMap mean stddev max min median mode skewness numberofvoxels volume");
	// Color code nuc by intensity mean ch4
	//run("Intensity Measurements 2D/3D", "input=C4 labels=NucObjectsMap mean stddev max min median mode skewness numberofvoxels volume");
	ColorCodeFromTable(Feature, bgsName+"-intensity-measurements", 1, 1);
}

// Color Code the nuc image
function ColorCodeFromTable(FeatureName,TableName, saveFlag, closeFlag)
{
	waitForUser("To color code nuc by "+FeatureName+" select the table  "+TableName+"  , and select  "+FeatureName+" column\nSet Min/Max as desired, but it is advised to select smaller Min value than automatically suggested ");
	selectWindow("NucObjectsMap");
	run("Z Project...", "projection=[Max Intensity]");
	run("Assign Measure to Label");
	run("Fire");
	run("Calibration Bar...", "location=[Lower Left] fill=White label=Black number=5 decimal=2 font=12 zoom=1 overlay");
	name = "ColorCodeObjects_"+TableName+"_"+FeatureName;
	rename(name);
	saveAs("Tiff", resFolder+origNameNoExt+"_"+name+".tif");
	if (saveFlag == 1) 
		saveTable(resFolder, origNameNoExt, TableName, closeFlag);
}

// save and close the table
//---------------------------------------------------------
function saveTable(resFolder, origNameNoExt, table_name, closeFlag)
{
	selectWindow(table_name);
	//saveAs("Results", resFolder+origNameNoExt+"_"+table_name+".csv");
	Table.save(resFolder+origNameNoExt+"_"+table_name+".csv"); 
	if (closeFlag == 1) run("Close");  // To close non-image window
}


// Setup Quantification Channels and Prepare Channel-images 
//---------------------------------------------------------
function SetupChannelImages(origIm)
{
	selectWindow(origIm);
	// Setup Quantification Channels
	getDimensions(width, height, numChannels, numSlices, numFrames);
	nQuantChannels = 0;
	for (n = 1; n <= numChannels; n++)
	{
		if ((n != NucDNAChannel) && (n != NucSegChannel))
		{
			if (nQuantChannels == 0) 
			{
				QuantChannels = newArray(1);
				QuantChannels[0] = n;
			}
			else 
				QuantChannels = Array.concat(QuantChannels,n);

			nQuantChannels++;
		}
	}
	print("Setting Quant Channel, numChannels=", numChannels, " Channel for Segmentation: ", NucSegChannel, "Channel for Normalization:", NucDNAChannel, "Channels for Quantification: ");
	Array.print(QuantChannels);

	// Duplicate image, and create composite image for visualization and quality control
	rename("Raw");
	getVoxelSize(pixelWidth, pixelHeight, pixelDepth, unit);
	run("Duplicate...", "title=composite duplicate");
	run("Make Composite");
	run("Z Project...", "projection=[Max Intensity]");
	if (SaveMaxProject == 1)
		saveAs("Tiff", resFolder+origNameNoExt+"_MaxProject.tif");
	rename("MaxProject");
	
	//run("Channels Tool...");
	selectWindow("composite");
	Stack.setDisplayMode("color");
	for (n = 0; n < numChannels; n++)
	{
		Stack.setChannel(n);
		run(ChColor[n]);
	}
	Stack.setDisplayMode("composite");
	
	// Split channles, and use simple names 
	selectWindow("Raw");
	run("Split Channels");
	for (n = 1; n <= numChannels; n++)
	{
		selectWindow("C"+n+"-Raw");
		rename("C"+n);
	}

	// Background subtraction - 2D rolling ball 
	RollingBallRadiusPixels = round(RollingBallRadius / pixelWidth); 
	if (NucDNAChannel != -1) 
		BgsChannels = Array.concat(NucDNAChannel, QuantChannels);
	else
		BgsChannels = QuantChannels;
		
	print("applying Rolling ball with radius ",RollingBallRadiusPixels, " pixels to channels:");
	Array.print(BgsChannels);
	//for (n = 0; n <= nQuantChannels; n++)
	for (n = 1; n <= numChannels; n++)
	{
		selectWindow("C"+n);
		bgsName = "C"+n+"_Bgs";
		run("Duplicate...", "title="+bgsName+" duplicate");
		run("Subtract Background...", "rolling="+RollingBallRadiusPixels+" stack");	
	}
} // end of SetupChannelImages

//--------------------------------------------------------------------------------
	
/* Nuc Segmentation: 
 * rolling ball (applied earlier), blur, threshold (2D, stack histogram), 3D connected componnet analysis, add to 3D Roi Manager 
 * output is "NucObjectsMap" 
 */
function NucSegmentation()
{
	//selectWindow("C1");
	//selectWindow("C"+NucSegChannel);
	selectWindow("C"+NucSegChannel+"_Bgs");
	getVoxelSize(pixelWidth, pixelHeight, pixelDepth, unit);
	run("Duplicate...", "title=Nuc duplicate");
	//run("Gaussian Blur...", "sigma=2 stack");
	run("Gaussian Blur...", "sigma="+NucBlurSig+" stack");
	
	//setAutoThreshold("Otsu dark stack");
	//setAutoThreshold("MaxEntropy dark stack");
	setAutoThreshold(NucAutoThMethod+" dark stack");
	if (CheckThresholdValidity == 1)
		waitForUser("Check Threshold Validity");
	selectWindow("Nuc");
	getThreshold(lower, upper);
	print("Threshold:",lower, upper);
	// Apply the selected threshold to all the stack to generate a binary image
	//run("Convert to Mask", "method=Otsu background=Dark black");
	run("Convert to Mask", "method="+NucAutoThMethod+" background=Dark black");
	// Convert from um to Pixel
	MinNucSizeVoxels = MinNucSize / (pixelWidth * pixelHeight * pixelDepth); 
	MaxNucsizeVoxels = MaxNucSize / (pixelWidth * pixelHeight * pixelDepth); 
	if (ExcludeNucOnEdge == 1)
		run("3D Objects Counter", "threshold=128 slice=11 min.="+MinNucSizeVoxels+" max.="+MaxNucsizeVoxels+" exclude_objects_on_edges objects");
	else 
		run("3D Objects Counter", "threshold=128 slice=11 min.="+MinNucSizeVoxels+" max.="+MaxNucsizeVoxels+" objects");
	selectWindow("Objects map of Nuc");
	rename("NucObjectsMap");
	run("glasbey inverted");
	
	//run("3D Manager Options", "volume compactness fit_ellipse integrated_density mean_grey_value std_dev_grey_value feret minimum_grey_value maximum_grey_value centroid_(unit) bounding_box sync distance_between_centers=10 distance_max_contact=1.80 drawing=Contour");
	//run("3D Manager Options", "volume compactness fit_ellipse integrated_density mean_grey_value std_dev_grey_value feret minimum_grey_value maximum_grey_value centroid_(unit) bounding_box sync distance_between_centers=10 distance_max_contact=1.80 drawing=Contour use_0");	
	run("3D Manager Options", "volume compactness fit_ellipse integrated_density mean_grey_value std_dev_grey_value feret minimum_grey_value maximum_grey_value centroid_(unit) bounding_box sync distance_between_centers=10 distance_max_contact=1.80 drawing=Contour use_0 use_1");
	run("3D Manager");
	Ext.Manager3D_Reset();
	selectWindow("NucObjectsMap");
	Ext.Manager3D_AddImage();
} // end of NucSegmentation

//--------------------------------------------------------------------------------

// Saves Run Parameters 
function SavePrms(resFolder, origNameNoExt)
{
	PrmFile = resFolder+origNameNoExt+"_Measure3DNucPrms.txt";
	File.saveString("NucSegChannel="+NucSegChannel, PrmFile);
	File.append("", PrmFile); 
	File.append("NucDNAChannel="+NucDNAChannel, PrmFile); 
	File.append("NucBlurSig="+NucBlurSig, PrmFile); 
	File.append("RollingBallRadius="+RollingBallRadius, PrmFile); 
	File.append("MinNucSize="+MinNucSize, PrmFile); 
	File.append("MaxNucSize="+MaxNucSize, PrmFile); 
	File.append("NucAutoThMethod="+NucAutoThMethod, PrmFile); 
	File.append("ExcludeNucOnEdge="+ExcludeNucOnEdge, PrmFile); 
} // SavePrms

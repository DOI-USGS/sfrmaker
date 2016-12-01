####sfr_classes depedencies


* #####ESRI Arcpy
which **must be added to the path of your main python distribution**.

	**Adding Arcpy to the python path**:  	
	1) **Make a file called: Desktop10.pth**, with the following lines:  

	```
	C:\ArcGIS\Desktop10.2\arcpy  
	C:\ArcGIS\Desktop10.2\bin64  
	C:\ArcGIS\Desktop10.2\ArcToolbox\Scripts
	```
Notes:  
These lines tell python where to find arcpy and its associated packages/libraries. The second line may be "bin" for 32 bit or "bin64" for 64 bit. If you are using ArcMap 10.0 or 10.1, "Desktop10.2" in the above path needs to be modified accordingly.

	2) **Place this file where python can find it.** For Anaconda on Windows 7, this path should work (replacing ```aleaf``` with your username):

		C:\Users\aleaf\AppData\Local\Continuum\Anaconda
		
	The **Lib** subfolder in this folder or the **site-packages** folder within Lib may also work. Python checks these folders by default when looking for modules; within these folders, it looks for files with the extension **.pth**, and checks for additional paths within those files (i.e., the arcpy paths listed above).  

* ####ArcGIS 10.3:
 If you are using ArcGIS 10.3, you may need to copy this file (or its equivilant, depending on where the ArcPy-related python distribution is located):  
 
	```  
	C:\Python27\ArcGISx6410.3\Lib\site-packages\DTBGGP64.pth
	```
	to 
	
	```
	C:\Users\USERNAME\Anaconda 	
	```
	(or wherever your main python distribution is installed)

####Running sfr_classes
With the XML input file, the SFR_main.py script (calling the methods in sfr_classes) will produce two tables, **Mat1** and **Mat2**, with SFR reach and segment information.
  
1) **Setup XML input file** (see EXAMPLE.XML in \<InputFiles> section) to point to the above input datasets  
  
* check settings in \<GlobalSettings> section; make sure that \<preproc> is set to **True**  
* Set \<elevflag> to NHDPlus (this will set segment end elevations from NHDPlus COMID min/max elevations in the **elevslope** table)

2) Make sure that the "infile" variable in **SFR_main.py** (see SFR_main_EXAMPLE.py) points to the XML input file. Also make sure that the calls to classes are entered correctly.

3) **Run SFR_main.py** by typing *python SFR_main.py* at the command prompt  

4) If a "manual intervention" message is encountered in screen output,  

* **check the following files:  **
	* **fix_com_IDs.txt:** Lists stream segments that have multiple ends and starts; usually these are 		streams that are broken into mulitple parts by the grid boundary. 
	* **boundary_manual_fix_issues.txt:** Lists stream segments that don't have any connections to other 		segments.  
* **edit the offending segments** (COMIDs) in the shapefile specified under \<IntermediateFiles>\<intersect> in the XML input file (usually this means deleting the parts of multi-part COMIDs that are isolated by the grid boundary, and possibly deleting any unconnected segments).  
  
5) set \<preproc> in the XML input file to False (meaning the existing \<intersect> shapefile will be read in lieu of the preprocessing operations). Then **rerun SFR_main.py**.

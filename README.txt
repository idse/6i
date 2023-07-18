Spatial single cell analysis of proteins in 2D human gastruloids using iterative immunofluorescence
Current Protocols

This code can be downloaded without data from :
https://github.com/idse/6i
and with data from:
https://www.dropbox.com/sh/edz56lj4teyb1pq/AACPhzhPmPCoozVprSSsDttWa?dl=0

This code has been tested on the following systems:

Matlab
	Requires Image Processing Toolbox, Computer Vision Toolbox
	Tested on:
	2021a
		- MacOS 10.15.7
	2022a 
		- MacOS 12.6.3
	2023a
		- Windows 11

Python v3.8.0
For Macs:
	Mplex_analysis environment
		- MacOS 10.15.7
		- MacOS 12.6.3
For Windows:
	Mplex_analysis_windows
		- Windows 10
	cellpose_environment
		- Windows 10

How to create python environments

Macs: Mplex_analysis environment installation
	Install anaconda
	open terminal 
	Navigate to “code” folder	cd /[Folder]/code
	enter: conda env create -f Mplex_analysis.yml
	(Cellpose included in Mplex_analysis.yml)
Windows: Mplex_analysis environment installation 
	Install anaconda
	open anaconda prompt
	Navigate to “code” folder	cd /[Folder]/code
	enter: conda env create -f Mplex_analysis_windows.yml
	enter: conda env create -f cellpose_environment.yml
	
Image stitching:
For your own experiment, copy stitchImages.m into each round subdirectory, run this script from each subdirectory. Script assumes 6x4 grid pattern

Alignment:
	1.	Run AlignRounds.m
	2.	Optional: looks at overlays with makeIFOverlays.m

Segmentation:
ilastik - (1.3.3post3-OSX, or 1.4.0 OSX) see instructions in paper
Cellpose - copy segmentation_cellpose.ipynb into RD1 data folder
		** for Windows create cellpose environment separately 

	1.	Execute ilastik segementation
	2.	Execute cellpsoe segmentation (for Windows activate cellpose environment, for Mac activate Mplex_analysis)
	3.	Run combineSegmentations.m

Quantification
	1.	Run SingleCellQuantification.m
	⁃	script will output cdv into data subfolder

Single Cell Analysis

Activate environment
	enter: conda activate Mplex_analysis
	(enter: conda deactivate to deactivate environment)

Open  MultiplexedSCAnalysis.ipynb, run code in order

Content of code folder

** asterisks indicate instruction to copy files when running tutorial on tutorial images

Code                              folder
|
|— AlignRounds.m	
|— cellpose_environment.yml
|— combineSegmentations.m
|— makeIFOverlays.m	
|— Mplex_analysis.yml
|— Mplex_analysis_windows.yml
|— MultiplexedSCAnalysis.ipynb
|— segmentation_cellpose.ipynb
|— singleCellQuantification.m
|— stitchImages.m
|— Transfer.m
|-- classes                       folder
|   |-- Colony.m
|   |-- Metadata.m
|   |-- Position.m
|-- data                             folder
|   |-- 230110_6i9rd_exp20_RD1_SMAD23_pAKT_SOX17               folder
|         |-- [230110_6i9rd_exp20_RD1_SMAD23_pAKT_SOX17.nd2] x 24      images
|         |-- ** COPY stitchImages.m here
|         |-- ** COPY segmentation_cellpose.ipynb
|   |-- 230111_6i9rd_exp20_RD2_GATA3_OTX2_LEF1 folder
|         |-- [230111_6i9rd_exp20_RD2_GATA3_OTX2_LEF1.nd2] x 24             images
|         |-- ** COPY stitchImages.m here
|   |-- sample_segmentations folder
|         |-- [images]
|   |-- 230713_CP_alldata.csv   
|-- functions                           folder
|   |-- external                        folder
|   |-- preprocessing              folder
|   |-- segmentation               folder
|   |-- cleanSubplot.m
|   |-- figurePosition.m
|   |-- readIntensityValues.m
|   |-- renameDuplicateChannels.m
|   |-- seglim.m
|   |-- stitchedlim.m
|-- processed_data                           folder
|   |-- 230110_6i9rd_exp20_RD1_SMAD23_pAKT_SOX17  folder
|        |—MIP                                                                          folder,  contains images
|   |--230111_6i9rd_exp20_RD2_GATA3_OTX2_LEF1         folder
|        |—MIP                                                                          folder,  contains images
|   |--230115_6i9rd_exp20_RD6_BCAT_TBX6_SMAD2       folder
|        |—MIP                                                                          folder, contains images



# HMM Analysis on Wakeman & Henson Dataset

### Prerequisites

To run these analyses you will require:

	- A UNIX-type computer system
	- FSL version 5.0.9 (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)
	- MatLab 2014b or greater
	- SPM12 MatLab Toolbox (http://www.fil.ion.ucl.ac.uk/spm/software/download/)
	- OSL2 MatLab Toolbox (https://ohba-analysis.github.io/osl-docs/)
	- HMM-MAR MatLab toolbox (https://github.com/OHBA-analysis/HMM-MAR)

OSL2 and HMM-MAR are provided as part of this download. SPM12 should be installed into the osl2 directory.

### Getting Started

The download containing core scripts and associated toolboxes can be found on our OSF page here: https://osf.io/ugjbr/

Please download the HMM\_Task\_Download.zip file. It is recommended to make sure the scripts subdirectory in the download folder contains the latest versions of the scripts from this repository.

Note: All folder paths referred to here are relative to the location of this download, ie all paths assume that your current directory is the top-level of this download.

Once the download is complete take a look into the toolboxes folder of this download. You will find a copy of HMM-MAR (github repo: https://github.com/OHBA-analysis/HMM-MAR) and OSL (docs: https://ohba-analysis.github.io/osl-docs/). These are the core MatLab toolboxes that will run the analyses.

Next, take a look into the scripts/ folder within the HMM download. This contains the main analysis script and a utils subfolder containing a small helper module with some useful functions used within the analysis.

Get the data:

1) Download revision 0.1.1 of the Wakeman and Henson dataset (http://doi.org/10.1038/sdata.2015.1). These files can be downloaded anywhere on your system which is convenient and that can be accessed through your MatLab file browser.

Install SPM12 and OSL (note these steps can be carried out whilst the data is downloading):

2) Download SPM12 and install into the toolboxes/osl folder, SPM is available here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/

3) Install OSL. Additional information about the dependencies and install path for OSL can be found in toolboxes/osl/README.md and toolboxes/osl/osl-core/README.md

Configure the analysis:

4) Open scripts/+utils/get\_studyinfo.m. This script defines the locations of the downloaded data, analysis scripts and the location to store results and is used throughout the analysis to load and save data. Edit the file paths to match the chosen locations on your system.

5) open scripts/hmm\_0\_initialise.m . This script adds all the relevant paths to your MatLab session and loads the information defined within get\_studyinfo.m . Edit the download\_path at the top to the location of this download.

Check installation

6) Open MatLab 2014b (or greater) and run hmm\_0\_initialise.m toolboxes to the path and load the information from get\_studyinfo.m .

7) Run utils.check\_raw\_data in your MatLab session. This should confirm that 'All raw data is found in studydir'. If not, check that the studydir path in scripts/+utils/get\_studyinfo.m matches the location of your downloaded data.

8) Run osl\_check\_installation in your MatLab session. This will print a range of information about your OSL and MatLab environment and highlight any missing dependencies.

### Main Analysis

The HMM analysis is completed in a MatLab session in which the hmm\_0\_initialise.m script has been excecuted.

The HMM analysis proceeds by following the hmm\_\* scripts within the scripts/ subfolder of this download. They should be run in order, though each script will define and load everything it needs. So you can exit matlab and start a new session in between running two of the scripts as long as hmm\_0\_initialise.m is excecuted in each new session.

The hmm\_\* scripts contain additional information and explanations which help to describe the preprocessing and analysis as it is run.

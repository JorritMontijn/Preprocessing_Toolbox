%runPostDetectionMaster Runs data post-preprocessing
%
%This preprocessing toolbox was created to minimize the effort of
%preprocessing two-photon calcium data and was created by Jorrit Montijn at
%the University of Amsterdam
%
%To finish all pre-processing steps, you should run the following scripts
%in order:
%
%(1) runPreProMaster
%	Use this function to pre-preprocess data; this includes meta data
%	retrieval (doPrePro), image processing (doImagePrePro) and pixel-based
%	stimulus response maps (doCalcPixelResponsiveness)
%(2) runDetectCells
%	This function is used to select regions of interest based on the
%	average images generated by the previous step; it includes broad
%	functionality including neuronal subtype differentation based on 960nm
%	reference images, pixel-based responsiveness maps, 
%	retrieval (doPrePro), image processing (doImagePrePro) and pixel-based
%	stimulus response maps (doCalcPixelResponsiveness), automatic border
%	detection for OGB and GCaMP, and across-recording alignment functions,
%	such as the custom-built automatic recursive locally affine subfield
%	reregistration algorithm (by JM), ROI shift detection, etc.
%(3) runPostDetectionMaster
%	This function runs all post-preprocessing steps consecutively for the
%	defined session/recordings; this includes the following steps:
%	runImageToTimeseries to extract fluorescence time traces from the image
%	data for all ROIs defined in the previous step;
%	runDetectCalciumTransients to calculate dF/F0 values from these traces
%	and detect spikes from the calcium transient data; and
%	runBuildSesFromPrePro to transform the preprocessing data format to a
%	more useable "ses" data format on which all data processing functions
%	are based
%(4) Optional: eye-tracking video analysis and Z-drift checking
%	(runDriftChecker/runAcrossRecordingDriftChecker) to confirm z-stability
%	of your recordings 
%
%	Version history:
%	1.0 - September 14 2012
%	Created by Jorrit Montijn
%	2.0 - May 16 2014
%	Updated help & comments [by JM]

clear all;

%set recordings
for i=1
	clearvars -except i
	if i == 1
		strSession = '20140423';
		vecRecordings = 1:2;
	end
	
	runImageToTimeseries;
	
	runDetectCalciumTransients;
	
	runBuildSesFromPrePro;
end
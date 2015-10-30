%runDriftChecker Runs z-drift checking algorithm
%Use this function to check the stability of your recordings over time. It
%requires a xyz stack as well as a time series recordings; it then
%calculates the similarity for each frame to the different planes in the
%z-stack, allowing you to investigate changes in z-level over time within a
%recording. It uses the subfunction doDriftCheck to perform the
%registration and the subfunction doPlotDriftCheck to show the results.
%doDriftCheck also saves a data file to the hard drive in the recording
%folder for later subsequent analysis
%
%	Version history:
%	1.0 - May 16 2014
%	Created by Jorrit Montijn


%% source data
strMasterPath = 'D:\Data\Processed\imagingdata\';
boolDoPlots = false;
	
%% input variables
for intMouse=[7 8] %[3 4 5 7 8]
	if intMouse == 3
		strDirRawStackZ = 'G:\Data\Raw\imagingdata\20140425\xyz_2';
		strSession = '20140425';
		vecRecordings = 1:8;
		vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
	elseif intMouse == 4
		strDirRawStackZ = 'G:\Data\Raw\imagingdata\20140507\xyz';
		strSession = '20140507';
		vecRecordings = 1:3;
		vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
	elseif intMouse == 5
		strDirRawStackZ = 'G:\Data\Raw\imagingdata\20140530\xyz';
		strSession = '20140530';
		vecRecordings = 1:7;
		vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
	elseif intMouse == 7
		strDirRawStackZ = 'G:\Data\Raw\imagingdata\20140711\xyz';
		strSession = '20140711';
		vecRecordings = 1:8;
		vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
	elseif intMouse == 8
		strDirRawStackZ = 'G:\Data\Raw\imagingdata\20140715\xyz';
		strSession = '20140715';
		vecRecordings = 1:4;
		vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
	end
	
	for intRecording=1:length(vecRecordings)
		%generate recording name
		strRecording = sprintf('xyt%02d',intRecording);
		
		%put data in structure
		sSourceData.strMasterPath = strMasterPath;
		sSourceData.strDirRawStackZ = strDirRawStackZ;
		sSourceData.strSession = strSession;
		sSourceData.strRecording = strRecording;
		
		%do drift check
		sRegZ = doDriftCheck(sSourceData);
		
		%plot output
		if boolDoPlots,doPlotDriftCheck(sRegZ);end
	end
end
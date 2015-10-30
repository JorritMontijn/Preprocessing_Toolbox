%runAcrossRecordingDriftChecker Runs z-drift checking algorithm over recordings
%Use this function to check the z-stability over your recordings. It
%requires a xyz stack as well as multiple time series recordings; it then
%calculates the similarity for at the start, middle and end of each
%recording to the different planes in the z-stack, allowing you to
%investigate slow drifts in z-level over time between recordings. It uses
%the subfunction doAcrossRecordingDriftCheck to perform the registration
%and the subfunction doPlotAcrossRecordingDriftCheck to show the results.
%doAcrossRecordingDriftCheck also saves a data file to the hard drive in
%the recording folder for later subsequent analysis
%
%	Version history:
%	1.0 - May 16 2014
%	Created by Jorrit Montijn

%% input variables
%source data
clear all
strMasterPath = 'D:\Data\Processed\imagingdata\';
boolDoPlots = true;

for intMouse=[3 4 5 7 8]
	if intMouse == 3
		strDirRawStackZ = 'G:\Data\Raw\imagingdata\20140425\xyz_2';
		strSession = '20140425';
		vecRecordings = 1:8;
		vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
	elseif intMouse == -1
		strDirRawStackZ = 'G:\Data\Raw\imagingdata\20140430\xyz';
		strSession = '20140430';
		vecRecordings = 1:8; %xyt01 is not used for analysis
		vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
	elseif intMouse == 4
		strDirRawStackZ = 'G:\Data\Raw\imagingdata\20140507\xyz';
		strSession = '20140507';
		vecRecordings = 1:3;
		vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
	elseif intMouse == 5
		strDirRawStackZ = 'G:\Data\Raw\imagingdata\20140530\xyz';
		strSession = '20140530';
		vecRecordings = 2:8;
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
	
	%put data in structure
	sSourceData.strMasterPath = strMasterPath;
	sSourceData.strDirRawStackZ = strDirRawStackZ;
	sSourceData.strSession = strSession;
	sSourceData.vecRecordings = vecRecordings;
	
	%do drift check
	%check for data
	strFile = [strMasterPath strSession filesep strSession '_across_rec_zreg.mat'];
	if exist(strFile,'file')
		load(strFile);
		sRegZ = struct;
		sRegZ.matRegistrationZ = matRegistrationZ;
		sRegZ.sData = sData;
		sRegZ.strSession = strSession;
		sRegZ.vecRecordings = vecRecordings;
		sRegZ.sRec = sRec;
		fprintf('Loaded pre-processed data %s; time is [%s]\n',strFile,getTime)
	else
		sRegZ = doAcrossRecordingDriftCheck(sSourceData);
	end
	
	%plot output
	if boolDoPlots
		if ~exist('sRegZ','var')
			%make structure
			sRegZ = struct;
			sRegZ.matRegistrationZ = matRegistrationZ;
			sRegZ.sData = sData;
			sRegZ.strSession = strSession;
			sRegZ.vecRecordings = vecRecordings;
			sRegZ.sRec = sRec;
		end
		doPlotAcrossRecordingDriftCheck(sRegZ);
		
		%save figure
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		strPath = ['D:\Data\Results\stimdetection\' strSession filesep];
		strFig = sprintf('%sacross_rec_zreg%s_raw',strPath,strSession);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
		close all;
	end
end
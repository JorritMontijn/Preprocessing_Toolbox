clear all;
for intSes=7:8
	clear sRec;
	clear strSession;
	clear vecRecordings;
	clear cellName;
	clear cellRefPaths;
	clear cellStimLog;
	
	sMD = struct; %structMetaData
	sMD.strMasterDir = 'D:\Data';
	sMD.strImgSource = '\Raw\imagingdata\';
	sMD.strLogSource = '\Raw\imaginglogs\';
	sMD.strImgTarget = '\Processed\imagingdata\';
	sMD.strLogTarget = '\Processed\imaginglogs\';
	sMD.strTemp = '\Temp\';
	if intSes==1
		strSession = '20130612';
		vecRecordings = [1 2];
	elseif intSes==2
		strSession = '20130625';
		vecRecordings = [1 2];
	elseif intSes==3
		strSession = '20131016';
		vecRecordings = [1:6];
	elseif intSes==4
		strSession = '20131022';
		vecRecordings = [1:3];
	elseif intSes==5
		strSession = '20140129';
		vecRecordings = [1:5];
	elseif intSes==6
		strSession = '20140314';
		vecRecordings = [8 9];
	elseif intSes==7
		strSession = '20140423';
		vecRecordings = [1:3];
	elseif intSes==8
		strSession = '20140425';
		vecRecordings = [9 10];
		
		%stim detect
	elseif intSes == 11
		strSes = '20140207';
		vecBlock = 1:8; %define whether neurons are in same population or not
	elseif intSes == 12
		strSes = '20140314';
		vecBlock = [1 1 1 1 1 1 1]; %define whether neurons are in same population or not
	elseif intSes == 13
		strSes = '20140425';
		vecBlock = [1 1 1 1 1 1 1 1]; %define whether neurons are in same population or not
	elseif intSes == 14
		strSes = '20140507';
		vecBlock = [1 1 1]; %define whether neurons are in same population or not
	elseif intSes == 15
		strSes = '20140530';
		vecBlock = [1 1 1 1 1 1 1]; %define whether neurons are in same population or not
	elseif intSes == 16
		strSes = '20140604';
		vecBlock = [1 1 1 1]; %define whether neurons are in same population or not
	elseif intSes == 17
		strSes = '20140711';
		vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
	elseif intSes == 18
		strSes = '20140715';
		vecBlock = [1 1 1 1]; %define whether neurons are in same population or not
	end
	
	
	%get data
	% define general metadata
	sPS = loadDefaultSettingsPrePro();%structProcessingSettings
	strMasterDir = 'D:\UvA_Backup\Data';
	strTargetDir = '\Processed\imagingdata\';
	
	% create filenames
	for intRec=vecRecordings
		strMatRec{intRec} = [strMasterDir strTargetDir strSession filesep sPS.strRecording{intRec} filesep strSession sPS.strRecording{intRec} '_prepro.mat']; %timeseries
	end
	
	%loop through recordings
	for intRec=vecRecordings
		
		sRecLoad = load(strMatRec{intRec});
		sRec = sRecLoad.sRec;
		clear sRecLoad;
		%process
		strNewDir = [sMD.strMasterDir sMD.strImgTarget sRec.strSession filesep sRec.sProcLib.strRecording filesep];
		sRec.sMD.strMasterDir = strMasterDir;
		sRec = doCalcPixelResponsiveness(sRec,strNewDir);
	end
end
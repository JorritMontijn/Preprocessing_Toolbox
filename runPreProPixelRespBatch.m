clear all;
for intSes=12
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
		strSession = '20120718';
		vecRecordings = [1 3];
	elseif intSes==2
		strSession = '20120720';
		vecRecordings = [1 3 4];
	elseif intSes==3
		strSession = '20121207';
		vecRecordings = [1 2];
	elseif intSes==4
		strSession = '20121212';
		vecRecordings = [1 3];
	elseif intSes==5
		strSession = '20130307';
		vecRecordings = [1 4 5];
	elseif intSes==6
		strSession = '20130313';
		vecRecordings = [2 3 5];
	elseif intSes==7
		strSession = '20130315';
		vecRecordings = [1 3];
	elseif intSes==8
		strSession = '20130612';
		vecRecordings = [2];
	elseif intSes == 81
		strSession = '20130627';
		vecRecordings = 1:4;
	elseif intSes==9
		strSession = '20130625';
		vecRecordings = [2];
	elseif intSes==10
		strSession = '20131016';
		vecRecordings = [1 3 5 6];
	elseif intSes==11
		strSession = '20131022';
		vecRecordings = [1 3];
	elseif intSes==12
		strSession = '20140129';
		vecRecordings = [2 3 4 5];
	end
	
	%get data
	% define general metadata
	sPS = loadDefaultSettingsPrePro();%structProcessingSettings
	strMasterDir = 'D:\Data';
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
		sRec = doCalcPixelResponsiveness(sRec);
	end
end
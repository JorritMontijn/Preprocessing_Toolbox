% initialization
%clear all;

%% input
if ~exist('strSession','var')
	strSession = '20130627';
	vecRecordings = 1:2;
end

%% define general metadata
sPS = loadDefaultSettingsPrePro();%structProcessingSettings
strMasterDir = 'D:\Data';
strTargetDir = '\Processed\imagingdata\';
sPS.boolUseParallel = false;

%% init parallel
if sPS.boolUseParallel && matlabpool('size') == 0
	matlabpool('open',4);
end

%% create filenames
for intRec=vecRecordings
	strMatROI{intRec} = [strMasterDir strTargetDir strSession filesep sPS.strRecording{intRec} filesep strSession sPS.strRecording{intRec} '_CD.mat']; %masks
	strMatRec{intRec} = [strMasterDir strTargetDir strSession filesep sPS.strRecording{intRec} filesep strSession sPS.strRecording{intRec} '_CD.mat']; %timeseries
end

%% load data
for intRec=vecRecordings
	sROI = load(strMatROI{intRec});
	if isfield(sROI.sRec,'CDG')
		sDC = doTransformCDGtosDC(sROI.sRec.CDG); %backward compatibility
	else
		sDC = sROI.sRec.sDC;
	end
	clear sROI;
	
	sRecLoad = load(strMatRec{intRec});
	sRec = sRecLoad.sRec;
	%clear sRecLoad;
	if ~isfield(sRec,'sDC')
		sRec.sDC = sDC;
		sRec.sMD.strMasterDir= 'D:\Data';
		sRec.sMD.strImgSource= '\Raw\imagingdata\';
		sRec.sMD.strLogSource= '\Raw\stimulationlogs\imaging\';
		sRec.sMD.strImgTarget= '\Processed\imagingdata\';
		sRec.sMD.strLogTarget= '\Processed\stimulationlogs\';
		sRec.sMD.strTemp= '\Temp\';
	end
	
	% run
	[saveFileName,sRec] = doImageToTimeseries(sRec,sDC,strMatRec{intRec});
end


%% close parallel workers if necessary
if sPS.boolUseParallel && matlabpool('size') > 0
	matlabpool('close');
end
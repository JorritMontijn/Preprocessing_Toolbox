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
strSourceMaster = 'D:\Data';

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
		sRec.sMD.strMasterDir= 'C:\Data';
		sRec.sMD.strImgSource= '\Raw\imagingdata\';
		sRec.sMD.strLogSource= '\Raw\stimulationlogs\imaging\';
		sRec.sMD.strImgTarget= '\Processed\imagingdata\';
		sRec.sMD.strLogTarget= '\Processed\stimulationlogs\';
		sRec.sMD.strTemp= '\Temp\';
	end
	
	% run
	strImPath = [strSourceMaster sRec.sMD.strImgTarget sRec.strSession filesep sRec.sProcLib.strRecording filesep];
	[saveFileName,sRec] = doImageToTimeseries(sRec,sDC,strMatRec{intRec},strImPath);
end


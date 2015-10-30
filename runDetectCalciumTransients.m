%runDetectCalciumTransients
%clear all;

%% input
if ~exist('strSession','var')
	strSession = '20140314';
	vecRecordings = 1:9;
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
strMasterPath = [strMasterDir strTargetDir strSession filesep];
strOldPath = cd(strMasterPath);
for intRec=vecRecordings
	strSubPath{intRec} = [sPS.strRecording{intRec} filesep]; %masks
end

%% run
for intRec=vecRecordings
	%load file
	strPath = [strMasterPath strSubPath{intRec}];
	cd(strPath);
	strMatRec = sprintf('%sxyt%02d_TS.mat',strSession,intRec);
	sLoad = load([strPath strMatRec]);
	sRec = sLoad.sRec;
	clear sLoad;
	
	%assign output
	sRec.spikeData = doDetectCalciumTransients(sRec);
	
	%get name&save
	strMatRec = sprintf('%sxyt%02d_TS.mat',strSession,intRec);
	if strcmp(getFlankedBy(strMatRec,'_','.mat'),'SD')
		%filename already has spikedetection append
		saveFileName = strMatRec(1:end-4);
	elseif strcmp(getFlankedBy(strMatRec,'_','.mat'),'TS')
		%filename ends with _TS.mat
		saveFileName = [strMatRec(1:end-7) '_SD.mat'];
	else
		saveFileName = [strMatRec(1:end-4) '_SD.mat'];
	end
	save(saveFileName,'sRec')
end
cd(strOldPath);

%% close parallel workers if necessary
if sPS.boolUseParallel && matlabpool('size') > 0
	matlabpool('close');
end
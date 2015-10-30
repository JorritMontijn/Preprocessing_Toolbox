%runBuildSesFromPrePro
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

%% create filenames
strMasterPath = [strMasterDir strTargetDir strSession filesep];
strOldPath = cd(strMasterPath);
for intRec=vecRecordings
	strRec{intRec} = [sPS.strRecording{intRec}]; %masks
end

%% run
for intRec=vecRecordings
	%get path
	strPath = [strMasterPath strRec{intRec} filesep];
	cd(strPath);
	strMatRec = sprintf('%s%s_SD.mat',strSession,strRec{intRec});


	%perform conversion
	ses = buildSesFromPrePro(strMatRec);

	%save
	if strcmp(getFlankedBy(strMatRec,'_','.mat'),'ses')
		%filename already has timeseries append
		saveFileName = strMatRec(1:end-4);
	elseif strcmp(getFlankedBy(strMatRec,'_','.mat'),'SD')
		%filename ends with _SD.mat
		saveFileName = [strMatRec(1:end-7) '_ses.mat'];
	else
		saveFileName = [strMatRec(1:end-4) '_ses.mat'];
	end
	save(saveFileName,'ses')
end
cd(strOldPath);
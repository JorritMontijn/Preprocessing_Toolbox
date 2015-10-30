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
for intMouse=[5 7 8]
	if intMouse == 3
		strSession = '20140425';
		vecRecordings = 1:8;
	elseif intMouse == 4
		strSession = '20140507';
		vecRecordings = 1:3;
	elseif intMouse == 5
		strSession = '20140530';
		vecRecordings = 1:7;
	elseif intMouse == 7
		strSession = '20140711';
		vecRecordings = 1:8;
	elseif intMouse == 8
		strSession = '20140715';
		vecRecordings = 1:4;
	end
	
	%create empty array
	vecDriftZ = [];
	for intRecording=1:length(vecRecordings)
		%msg
		fprintf('Starting %s%s [%s]\n',strSession,sprintf('xyt%02d',intRecording),getTime);
		
		%generate recording name
		strRecording = sprintf('xyt%02d',intRecording);
		strRecPath = [strMasterPath strSession filesep strRecording filesep];
		strFilename = sprintf('%s%s_zreg.mat',strSession,strRecording);
		
		%load data
		sRegZ = load([strRecPath strFilename]);
		
		if isfield(sRegZ,'sRegZ')
			sRegZ = sRegZ.sRegZ;
			vecThisDrift = sRegZ.vecDriftZ;
		else
			%get z-shifts
			close;
			vecThisDrift = doPlotWithinRecordingDriftCheck(sRegZ,true);
			sRegZ.vecDriftZ = vecThisDrift;
			
			%save data & figure
			save([strRecPath strFilename],'sRegZ');
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			export_fig([strRecPath strSession strRecording '_DriftZ.tif']);
			export_fig([strRecPath strSession strRecording '_DriftZ.pdf']);
		end
		vecDriftZ = [vecDriftZ vecThisDrift]; %#ok<AGROW>
	end
	%save data
	save([strMasterPath strSession filesep  strSession '_zdrift_aggregate.mat'],'vecDriftZ');
end
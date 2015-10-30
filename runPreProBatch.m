clear all;
for intSes=1:7
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
		
		%assign filenames
		cellName{1} = sprintf('OT1');
		cellRefPaths{1} =[sMD.strMasterDir sMD.strImgSource strSession filesep cellName{1} filesep];
		cellName{3} = sprintf('OT2');
		cellRefPaths{3} =[sMD.strMasterDir sMD.strImgSource strSession filesep cellName{3} filesep];
		
		cellStimLog{1} = '20120718_MP_PilotIII_OrientationTuning_OT1.mat'; %name of the stimulation log
		cellStimLog{3} = '20120718_MP_PilotIII_OrientationTuning_OT2.mat';
	elseif intSes==2
		
		strSession = '20120720';
		vecRecordings = [4];
		
		cellName{1} = sprintf('xyt01_OT1');
		cellRefPaths{1} =[sMD.strMasterDir sMD.strImgSource strSession filesep cellName{1} filesep];
		cellName{3} = sprintf('xyt03_OT2');
		cellRefPaths{3} =[sMD.strMasterDir sMD.strImgSource strSession filesep cellName{3} filesep];
		cellName{4} = sprintf('xyt04_OT3');
		cellRefPaths{4} =[sMD.strMasterDir sMD.strImgSource strSession filesep cellName{4} filesep];
		
		cellStimLog{1} = '20120720_MP_PilotIII_OrientationTuning_24_OT1.mat'; %name of the stimulation log
		cellStimLog{3} = '20120720_MP_PilotIII_OrientationTuning_OT2.mat';
		cellStimLog{4} = '20120720_MP_PilotIII_OrientationTuning_OT3.mat';
		
	elseif intSes==3
		strSession = '20121207';
		vecRecordings = [1 2];
		
		%assign filenames
		for intRec=vecRecordings
			cellName{intRec} = sprintf('xyt%02d',intRec);
			cellRefPaths{intRec} =[sMD.strMasterDir sMD.strImgSource strSession filesep cellName{intRec} filesep];
		end
		cellStimLog{1} = '20121207xyt01_MP_PilotIII_OrientationTuning_28.mat';
		cellStimLog{2} = '20121207xyt02_MP_PilotIII_OrientationTuning_28_50Hz.mat';
	elseif intSes==4
		strSession = '20121212';
		vecRecordings = [1 3];
		
		%assign filenames
		for intRec=vecRecordings
			cellName{intRec} = sprintf('xyt%02d',intRec);
			cellRefPaths{intRec} =[sMD.strMasterDir sMD.strImgSource strSession filesep cellName{intRec} filesep];
		end
		cellStimLog{1} = '20121212xyt01_MP_PilotIII_OrientationTuning_30.mat';
		cellStimLog{3} = '20121212xyt03_MP_PilotIII_OrientationTuning_30_2.mat';
	elseif intSes==5
		strSession = '20130307';
		vecRecordings = [1 4 5];
		
		%assign filenames
		for intRec=vecRecordings
			cellName{intRec} = sprintf('xyt%02d',intRec);
			cellRefPaths{intRec} =[sMD.strMasterDir sMD.strImgSource strSession filesep cellName{intRec} filesep];
		end
		cellStimLog{1} = '20130307_MP_OrientationTuning_35_anesth.mat';
		cellStimLog{4} = '20130307_MP_OrientationTuning_35_awake.mat';
		cellStimLog{5} = '20130307_MP_OrientationTuning_35_awake2.mat';
	elseif intSes==6
		strSession = '20130313';
		vecRecordings = [2 3 5];
		
		%assign filenames
		for intRec=vecRecordings
			cellName{intRec} = sprintf('xyt%02d',intRec);
			cellRefPaths{intRec} =[sMD.strMasterDir sMD.strImgSource strSession filesep cellName{intRec} filesep];
		end
		cellStimLog{2} = '20130313_MP_OrientationTuning_36_anesth.mat';
		cellStimLog{3} = '20130313_MP_OrientationTuning_36_awake.mat';
		cellStimLog{5} = '20130313_MP_OrientationTuning_36_LA4.mat';
	elseif intSes==7
		strSession = '20130315';
		vecRecordings = [1 3];
		
		%assign filenames
		for intRec=vecRecordings
			cellName{intRec} = sprintf('xyt%02d',intRec);
			cellRefPaths{intRec} =[sMD.strMasterDir sMD.strImgSource strSession filesep cellName{intRec} filesep];
		end
		cellStimLog{1} = '20130315_MP_OrientationTuning_anesth.mat';
		cellStimLog{3} = '20130315_MP_OrientationTuning_awake.mat';
	end
	runPreProMaster
end
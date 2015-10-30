%runBuildSesFromPrePro

function ses = buildSesFromPrePro(strMatRec)
	
	%define struct
	ses = struct;
	
	% load fluorescence data
	fprintf('Starting transformation from sRec type to ses type... [%s]\n',strMatRec);
	sLoad = load(strMatRec);
	sRec = sLoad.sRec;
	clear sLoad;
	
	%% get general info
	ses.strImPath = [sRec.sMD.strMasterDir sRec.sMD.strImgTarget sRec.strSession filesep sRec.sProcLib.strRecording filesep];
	ses.intLengthT = length(num2str(sRec.sProcLib.t-1));
	ses.strTargetIm = ['t%0' num2str(ses.intLengthT) 'd_ch%02d.tif'];
	if isfield(sRec, 'strMouse')
		ses.mouse = sRec.strMouse;
	else
		ses.mouse = sRec.strSession;
	end
	ses.session = sRec.strSession;
	ses.recording = sRec.intRecording;
	if sRec.sProcLog.boolXMLFound
		ses.xml = sRec.xml.sData;
		ses.date = sRec.xml.sData.strStartTime;
	else
		ses.xml = [];
		ses.date = sRec.strSession;
	end
	if sRec.sProcLog.boolStimLogFound
		ses.logname = sRec.sRawLib.strStimLog;
	else
		[sRec,intFlag] = loadStimLog(sRec);
		if intFlag ~= 0
			ses.logname = sRec.sRawLib.strStimLog;
		else
			ses.logname = '[]';
		end
	end
	fprintf('Extracting info from recording %d of session %s (date: %s; log: %s)\n',ses.recording,ses.session,ses.date,ses.logname);
	
	% calculate time
	ses.time.exp = ses.date;
	ses.time.start = ses.date;
	if sRec.sProcLog.boolXMLFound
		strT = sRec.xml.sData.strActualImageSizeT;
		intLocM = strfind(strT,'m');
		intMins = str2double(strT(1:(intLocM-1)));
		ses.time.dur = str2double(strT((intLocM+1):(end-1))) + intMins*60;
		ses.samplingFreq = sRec.xml.sData.dblFrameDur;
	else
		ses.time.dur = 0;
	end
	
	% calculate depth
	ses.depth = 0;
	
	% copy dimension sizes
	ses.size = sRec.sProcLib;
	
	% get sampling frequency
	frameTime = ses.time.dur/ses.size.t;
	ses.samplingFreq = 1/frameTime;
	
	% get anesthesia level
	ses.anesthesia = [];
	
	% get experiment description
	ses.experiment = ses.logname;
	ses.protocol = ses.logname;
	
	%% get stimuli
	if sRec.sProcLog.boolStimLogFound
		StimuliLoaded = 1;
		ses.structStim = getStimFramesPrePro(sRec);
	else
		StimuliLoaded = 0;
	end
	
	%% get movement correction output
	if isfield(sRec.sProcLog,'registration')
		for intCh = 1:numel(sRec.sProcLog.registration.Ch)
			ses.movementcorrection.ch(intCh).performance =sRec.sProcLog.registration.Ch(intCh).performance;
			ses.movementcorrection.ch(intCh).phaseCorr = sRec.sProcLog.registration.Ch(intCh).phaseCorr;
			ses.movementcorrection.ch(intCh).xReg = sRec.sProcLog.registration.Ch(intCh).xReg;
			ses.movementcorrection.ch(intCh).yReg = sRec.sProcLog.registration.Ch(intCh).yReg;
		end
	end
	
	%% get all objects
	intObjects = length(sRec.timeseries.roi);
	intTypes = length(sRec.sDC.metaData.cellType);
	vecObjectCounter = zeros(1,intTypes);
	for intObject = 1:intObjects
		
		%get type
		intType = sRec.sDC.ROI(intObject).intType;
		strObjectType = sRec.sDC.metaData.cellType{intType};
		
		%increment counter
		intThisCounter = vecObjectCounter(intType) + 1;
		vecObjectCounter(intType) = intThisCounter;
		
		%get dynamic field name
		strField = strObjectType;
		
		%set data
		ses.(strField)(intThisCounter).id = intThisCounter;
		ses.(strField)(intThisCounter).oldId = intObject;
		ses.(strField)(intThisCounter).x = sRec.sDC.ROI(intObject).intCenterX;
		ses.(strField)(intThisCounter).y = sRec.sDC.ROI(intObject).intCenterY;
		ses.(strField)(intThisCounter).matPerimeter = sRec.sDC.ROI(intObject).matPerimeter;
		ses.(strField)(intThisCounter).matMask = sRec.sDC.ROI(intObject).matMask;
		ses.(strField)(intThisCounter).size = sum(sRec.sDC.ROI(intObject).matMask(:));
		
		%get new data
		if isfield(sRec.sDC.ROI(intObject),'intPresence')
			if isscalar(sRec.sDC.ROI(intObject).intPresence),intPresence=sRec.sDC.ROI(intObject).intPresence;else intPresence = 1;end
			ses.(strField)(intThisCounter).intPresence = intPresence;
			ses.(strField)(intThisCounter).strPresence = sRec.sDC.metaData.cellPresence{intPresence};
			if isscalar(sRec.sDC.ROI(intObject).intRespType),intRespType=sRec.sDC.ROI(intObject).intRespType;else intRespType = 1;end
			ses.(strField)(intThisCounter).intRespType = intRespType;
			ses.(strField)(intThisCounter).strRespType = sRec.sDC.metaData.cellRespType{intRespType};
		end
		
		% get raw fluorescence
		ses.(strField)(intThisCounter).F = sRec.timeseries.roi(intObject).F;
		ses.(strField)(intThisCounter).npF = sRec.timeseries.roi(intObject).npF;
		
		% info from spikeData
		if ismember(intType,sRec.sDC.metaData.vecNeurons)
			ses.(strField)(intThisCounter).dFoF = sRec.spikeData(intObject).dFoF;
			ses.(strField)(intThisCounter).expFit = sRec.spikeData(intObject).expFit;
			ses.(strField)(intThisCounter).apFrames = sRec.spikeData(intObject).apFrames;
			ses.(strField)(intThisCounter).apSpikes = sRec.spikeData(intObject).apSpikes;
			if isfield(sRec.spikeData(intObject),'vecSpikes') && ~isempty(sRec.spikeData(intObject).vecSpikes)
				ses.(strField)(intThisCounter).vecSpikes = sRec.spikeData(intObject).vecSpikes;
			end
		end
	end
	
	%% transform interneurons to neuron subtype
	%assign new additional field to all neurons
	intNeurons = numel(ses.neuron);
	for intNeuron=1:intNeurons
		ses.neuron(intNeuron).type = 'neuron';
	end
	
	%convert interneuron structures to neurons
	cellInterneurons{1} = 'PV';
	cellInterneurons{2} = 'SOM';
	cellInterneurons{3} = 'VIP';
	
	for intInterneuron=1:length(cellInterneurons)
		strType = cellInterneurons{intInterneuron};
		
		if isfield(ses,strType)
			intCells = numel(ses.(strType));
			for intCell=1:intCells
				intNeuron = intNeuron + 1;
				structCell = ses.(strType)(intCell);
				
				cellFields = fieldnames(structCell);
				for intField=1:length(cellFields)
					ses.neuron(intNeuron).(cellFields{intField}) = structCell.(cellFields{intField});
				end
				ses.neuron(intNeuron).type = strType;
			end
		end
	end
	
	%% check all types and build counter fields
	for intType = 1:length(sRec.sDC.metaData.cellType)
		%get name & number of objects
		strType = sRec.sDC.metaData.cellType{intType};
		strField = sprintf('int_%s',strType);
		if isfield(ses,strType)
			intObjects = numel(ses.(strType));
			
			%assign to field
			ses.(strField) = intObjects;
			
			%for backward compatibility
			if strcmp(strField,'neuron')
				ses.nNeurons = intObjects;
			end
		end
	end
	ses.cellType = sRec.sDC.metaData.cellType;
	
	%% transform spikes to conform to dFoF if spikes are present
	if ~isempty(ses.neuron(1).expFit) && ~isfield(ses.neuron(1),'vecSpikes')
		disp('Transforming spike data to conform to dF/F format..');
		ses = doTransformSpikes(ses);
	end
end



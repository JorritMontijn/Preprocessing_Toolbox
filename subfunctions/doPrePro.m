function sRec = doPrePro(sRec)
	%help
	
	
	
	%starting message
	fprintf('\nStarting preprocessing of session %s (recording %s; name=%s)\n',sRec.strSession,sRec.sRawLib.strRecording,sRec.sRawLib.strName);
	
	%make log structure
	sRec.sProcLog = struct; %processing log
	
	%check settings
	if sRec.sPS.boolSupervised
		fprintf('Default values are:\n');
		disp(sRec.sPS)
		strInput = input('Is this correct? (Y/N)\n   ','s');
		if ~strcmpi(strInput,'Y')
			error([mfilename ':WrongSetting'],'Aborted due to incorrect settings');
		end
	end
	sRec.sProcLog.boolDefaultFound = true;
	
	% create processed data structure
	sRec.sProcLib = struct;
	sRec.sProcLib.strRecording = sRec.sPS.strRecording{sRec.intRecording};
	
	%is temp data present?
	sRec.sProcLog.boolTempPresent = checkTempPrePro(sRec);
	if sRec.sProcLog.boolTempPresent && sRec.sPS.boolSupervised
		warning([mfilename ':TempDataPresent'],'Temp data is present... Perhaps pre-processing failed last time.\n');
		strInput = input('Do you wish to continue anyway? (Y/N)\n   ','s');
		if ~strcmpi(strInput,'Y')
			error([mfilename ':TempDataPresent'],'Aborted due to temp data presence');
		end
	end
	
	%make path variables
	strOldDir = [sRec.sMD.strSourceDir sRec.sMD.strImgSource sRec.strSession filesep sRec.sRawLib.strRecording filesep];
	strOldFile = [sRec.sRawLib.strName '_txxxx_ch0x.tif'];
	strNewDir = [sRec.sMD.strMasterDir sRec.sMD.strImgTarget sRec.strSession filesep sRec.sProcLib.strRecording filesep];
	strNewFile = 'txxxx_ch0x.tif';
	
	%display old and new location
	fprintf('Will transfer %s%s to %simages%s%s during preprocessing\n',strOldDir,strOldFile,strNewDir,filesep,strNewFile);
	
	%load xml data
	if sRec.sPS.boolDoXMLRead
		if exist([strOldDir sRec.sRawLib.strName '_Properties.xml'],'file')
			strXMLfile = [strOldDir sRec.sRawLib.strName '_Properties.xml'];
			sRec.sProcLog.boolXMLFound = true;
		elseif exist([strOldDir 'MetaData' filesep sRec.sRawLib.strName '_Properties.xml'],'file')
			strXMLfile = [strOldDir 'MetaData' filesep sRec.sRawLib.strName '_Properties.xml'];
			sRec.sProcLog.boolXMLFound = true;
		else
			sRec.sProcLog.boolXMLFound = false;
		end
		
		if sRec.sProcLog.boolXMLFound == true
			fprintf('Reading XML file %s\n',strXMLfile);
			sRec.xml = struct;
			sRec.xml.strFileLocation = strXMLfile;
			[sRec.xml.sData,intOutFlag] = loadXMLPrePro(strXMLfile);
			if intOutFlag == 1
				fprintf('XML file succesfully loaded!\n')
			elseif intOutFlag == 0
				warning([mfilename ':XMLReadError'],'Variable retrieval incomplete\n');
			else
				error([mfilename ':XMLReadError'],'XML file loading failed\n');
			end
		else
			error([mfilename ':XMLFileNotFound'],sprintf('XML file not found at %s\n',strOldDir)); %#ok<SPERR>
		end
	else
		sRec.sProcLog.boolXMLFound = nan;
	end
	
	%check stimulation log data
	if sRec.sPS.boolDoStimLogRead
		strStimLog = [sRec.sMD.strSourceDir sRec.sMD.strLogSource sRec.strSession filesep sRec.sRawLib.strStimLog];
		if exist(strStimLog,'file')
			fprintf('Stimulation log was found!\n')
			sRec.sProcLog.boolStimLogFound = true;
			sRec.sStim = load(strStimLog);
		else
			error([mfilename ':NoStimLog'],'Stimulation log %s was not present\n',strStimLog);
		end
	else
		sRec.sProcLog.boolStimLogFound = nan;
	end
	
	%test-load unprocessed data
	% - load and check image consistency
	% - - make universal reference library to images
	if sRec.sPS.boolDoTestLoad
		fprintf('Checking image file consistency...\n')
		sImCheck = doCheckImageFilesPrePro(sRec);
		sRec.sProcLog.boolTestLoaded = true;
		
		sRec.sProcLib.t = sImCheck.t;
		sRec.sProcLib.z = sImCheck.z;
		sRec.sProcLib.x = sImCheck.x;
		sRec.sProcLib.y = sImCheck.y;
		sRec.sProcLib.ch = sImCheck.ch;
		fprintf('\b Done!\n')
	else
		sRec.sProcLog.boolTestLoaded = nan;
	end
	
	%check if image data conforms to XML data
	if sRec.sProcLog.boolTestLoaded && sRec.sProcLog.boolXMLFound
		if	 (sRec.xml.sData.intImageSizeX == sRec.sProcLib.x && ...
				sRec.xml.sData.intImageSizeY == sRec.sProcLib.y && ...
				sRec.xml.sData.intImageSizeT == sRec.sProcLib.t && ...
				sRec.xml.sData.intImageChannels == sRec.sProcLib.ch)
			fprintf('Image information is consistent with XML data\n');
		else
			warning([mfilename ':DataInconsistent'],'Image information is not consistent with XML data!')
			fprintf('Dimensions inconsistent: \nx: %d[XML]; %d[ProcLib]\ny: %d[XML]; %d[ProcLib]\nt: %d[XML]; %d[ProcLib]\nCh: %d[XML]; %d[ProcLib]\n\n',...
				sRec.xml.sData.intImageSizeX,sRec.sProcLib.x,...
				sRec.xml.sData.intImageSizeY,sRec.sProcLib.y,...
				sRec.xml.sData.intImageSizeT,sRec.sProcLib.t,...
				sRec.xml.sData.intImageChannels,sRec.sProcLib.ch);
		end
		%get actual duration of experiment
		strT = sRec.xml.sData.strActualImageSizeT;
		intLocM = strfind(strT,'m');
		intMins = str2double(strT(1:(intLocM-1)));
		dblTotDurSecs = str2double(strT((intLocM+1):(end-1))) + intMins*60;
		
		%get calcium channel
		sRec.sProcLib.CaCh = sRec.sProcLib.ch;
		
		%add values to library
		sRec.sProcLib.dblTotDurSecs = dblTotDurSecs;
	end
	
	% - ask for phase correction
	if sRec.sPS.boolDoPhaseCorrect
		sRec.sProcLog.dblMustPhaseCorrect = doCheckPhasePrePro(sRec);
	else
		sRec.sProcLog.dblMustPhaseCorrect = nan;
	end
	
	sRec.sProcLog.boolRawLibPresent = true;
	sRec.sProcLog.boolProcLibPresent = true;
	
	%make directories if not present
	structDir = dir(strNewDir);
	if numel(structDir) == 0;
		mkdir(strNewDir);
	end
	
	%save structure to temp directory
	strFile = [strNewDir sRec.strSession sRec.sProcLib.strRecording '_prepro.mat'];
	save(strFile,'sRec');
	fprintf('Saved recording structure to %s\n',strFile)
	
	%run pre-processing
	% - smoothing settings
	% - save session data (x, y, z, t, ch dimensions)
	% - load primary registration image as the average of the first 100 imgs
	% - loop through images
	% - - load image
	% - - update raw average image
	% - - apply phase correction (and save some examples)
	% - - apply smoothing (and save some examples)
	% - - apply img registration (and save some examples)
	% - - update processed average image
	% - update flags for processing steps taken
	% save images
	% - raw average
	% - HQ raw average
	% - processed average
	% - HQ processed average
	% - overlay of processed averages
	% - HQ overlay of processed averages
	
	
	
	
end


function [saveFileName,sRec] = doImageToTimeseries(sRec,sDC,strMatRec)
	%retrieves timeseries from image ROIs
	
	%define image variables
	if isfield(sRec.sProcLib,'CaCh')%which channel has calcium data?
		intCaCh = sRec.sProcLib.CaCh;
	else
		intCaCh = 1;
	end
	if isfield(sRec.sMD,'strMasterDir')
		strImPath = [sRec.sMD.strMasterDir sRec.sMD.strImgTarget sRec.strSession filesep sRec.sProcLib.strRecording filesep];
	else
		strImPath = [sRec.sMD.strImgTarget sRec.strSession filesep sRec.sProcLib.strRecording filesep];
	end
	intLengthT = length(num2str(sRec.sProcLib.t-1));
	strTargetIm = ['t%0' num2str(intLengthT) 'd_ch%02d.tif'];
	
	%check if files exist
	sDir = dir(strImPath);
	if isempty(sDir)
		vecSep = strfind(strMatRec,filesep);
		strImPath = strMatRec(1:vecSep(end));
		sRec.sMD.strMasterDir = strMatRec(1:(vecSep(2)-1));
	end
	
	%% processing
	%get masks per cell and per neuropil area
	fprintf('Creating masks, please wait...\n');
	[MaskCell, MaskCellNeuropil] = getROImasksPrePro(sRec,sDC);
	
	%pre-allocate
	intObjects = length(MaskCell);
	for intObject=1:intObjects
		timeseries.roi(intObject).F = zeros(1,sRec.sProcLib.t);
		timeseries.roi(intObject).npF = zeros(1,sRec.sProcLib.t);
	end
	fracPrev = -1;
	
	% calculate cell-timelines
	fprintf('Extracting timeseries, please wait...\n');
	
	%% output
	for intT=1:sRec.sProcLib.t
		%read images
		strIm = [strImPath 'images' filesep sprintf(strTargetIm,intT-1,intCaCh-1)];
		image = imread(strIm);
		matIm = im2double(image);
		
		for intObject=1:intObjects
			% get average roi-fluorescence per channel
			timeseries.roi(intObject).F(intT) = mean(matIm(MaskCell{intObject}));
			

			% get average surrounding neuropil-fluorescence per channel
			timeseries.roi(intObject).npF(intT) = mean(matIm(MaskCellNeuropil{intObject}));
		end
		
		%send msg
		intFrac = intT/sRec.sProcLib.t;
		fracNow = round(intFrac * 100);
		if fracPrev ~= fracNow
			tStamp = fix(clock);
			strPlace = sprintf('Processing %s%s... Now at %d%% [%02d:%02d:%02d]', sRec.strSession, sRec.sProcLib.strRecording, fracNow, tStamp(4),tStamp(5),tStamp(6));
			disp(strPlace);
			fracPrev = fracNow;
		end
	end
	
	% save cellData into sRec
	sRec.timeseries = timeseries;
	sRec.masks.MaskCell = MaskCell;
	sRec.masks.MaskCellNeuropil = MaskCellNeuropil;
	
	if strcmp(getFlankedBy(strMatRec,'_','.mat'),'TS')
		%filename already has timeseries append
		saveFileName = strMatRec(1:end-4);
	elseif strcmp(getFlankedBy(strMatRec,'_','.mat'),'CD')
		%filename ends with _CD.mat
		saveFileName = [strMatRec(1:end-7) '_TS.mat'];
	elseif strcmp(getFlankedBy(strMatRec,'_','.mat'),'prepro')
		saveFileName = [strMatRec(1:end-11) '_TS.mat'];
	else
		saveFileName = [strMatRec(1:end-4) '_TS.mat'];
	end
	
	save(saveFileName,'sRec')
	
end


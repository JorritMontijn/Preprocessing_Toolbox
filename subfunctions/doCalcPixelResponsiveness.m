function sRec = doCalcPixelResponsiveness(sRec,strNewDir,boolCalcPixCorr)
	%doCalcPixelResponsiveness Calculates pixel-based metrics for cell detection.
	%   Syntax: sRec = doCalcPixelResponsiveness(sRec,strNewDir)
	%
	%   Input:
	%	- sRec, pre-processing data structure
	%	Optional:
	%	- strNewDir: will save output to this path
	%	Output:
	%	- sRec, with added field "sPixResp"
	
	%% prep
	%get switches
	if ~exist('boolCalcPixCorr','var')
		boolCalcPixCorr = true;
	end
	
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
	
	% calculate cell-timelines
	fprintf('Extracting timeseries for pixel responsiveness, please wait...\n');
	
	%get stimulus timing
	structStim = getStimFramesPrePro(sRec);
	
	%pre-allocate look-up table and other variables
	vecOriLookup = getOriListFromTrials(structStim.Orientation);
	intStimTypes = length(vecOriLookup);
	intStimNumber = length(structStim.FrameOn);
	intReps = intStimNumber/intStimTypes;
	
	intBaselineIndex = intStimTypes + 1;
	intLastType = nan;
	intTraceLength = max(max(abs(structStim.FrameOn-structStim.FrameOff)),structStim.FrameOn(1));
	matTrace = nan(sRec.sProcLib.y,sRec.sProcLib.x,intTraceLength);
	intTraceCounter = 0;
	vecRepCounter = zeros(1,intBaselineIndex);
	
	%pre-allocate response matrix
	matPixelStimResp = nan(sRec.sProcLib.y,sRec.sProcLib.x,intStimTypes,intReps);
	matPixelNormResp = nan(sRec.sProcLib.y,sRec.sProcLib.x,intStimTypes,intReps);
	matPixelBaseResp = nan(sRec.sProcLib.y,sRec.sProcLib.x,1,intStimNumber);
	
	%pre-allocate 30-second buffer for pixel-wise correlations
	intSubY = floor(sRec.sProcLib.y/2);
	intSubX = floor(sRec.sProcLib.x/2);
	intBufferCounter = 0;
	intTemporalAveraging = 5;
	intTempAvgCounter = 0;
	intBufferT = 2000;
	intSpatialWindow = 12;
	intSurrSize = intSpatialWindow*2+1;
	intCorrsCounter = 0;
	matPixRespBuffer = nan(intSubY,intSubX,intBufferT,'single');
	matPixRespAverageBuffer = zeros(intSubY,intSubX,'double');
	matPixelCorrMat =  zeros(intSubY,intSubX,intSurrSize,intSurrSize,'single');
	
	%% perform stim-resp extraction
	fracPrev = -1;
	for intT=1:sRec.sProcLib.t
		%read image
		strIm = [strImPath 'images' filesep sprintf(strTargetIm,intT-1,intCaCh-1)];
		image = imread(strIm);
		matIm = im2double(image);
		
		%get stimulus at this frame
		dblOriAtFrame = getStimAtFrame(structStim,intT);
		if dblOriAtFrame == 999 %baseline
			intStimIndex = intBaselineIndex;
		else
			intStimIndex = find(vecOriLookup==dblOriAtFrame);
		end
		
		%check for change
		if intStimIndex ~= intLastType
			%change
			if ~isnan(intLastType)
				%if it's not very first frame
				
				%get repetition and increment for this type
				vecRepCounter(intLastType) = vecRepCounter(intLastType) + 1;
				intRep = vecRepCounter(intLastType);
				
				%check if stim or baseline
				if intBaselineIndex ~= intLastType
					%get mean response during this trace and put into output
					matThisResp = nanmean(matTrace,3);
					matPixelStimResp(:,:,intLastType,intRep) = matThisResp;
					
					matLastBase = matPixelBaseResp(:,:,1,vecRepCounter(intBaselineIndex));
					matPixelNormResp(:,:,intLastType,intRep) = matThisResp - matLastBase;
					
					%check if it's last to switch to baseline
					if intBaselineIndex == intStimIndex && vecRepCounter(intBaselineIndex) == intStimNumber
						break;
					end
				else
					%put in matrix
					matThisResp = nanmean(matTrace,3);
					matPixelBaseResp(:,:,1,intRep) = matThisResp;
				end
			end
			%clear trace buffer
			matTrace = nan(size(matTrace));
			intTraceCounter = 0;
			intLastType = intStimIndex;
		end
		
		%add current frame to trace buffer
		intTraceCounter = intTraceCounter + 1;
		matTrace(:,:,intTraceCounter) = matIm;
		if boolCalcPixCorr
			%add current frame to averaging buffer
			matSubImY = (matIm(1:2:end,:) + matIm(2:2:end,:))/2;
			matSubImXY = (matSubImY(:,1:2:end) + matSubImY(:,2:2:end))/2;
			intTempAvgCounter = intTempAvgCounter + 1;
			matPixRespAverageBuffer = matPixRespAverageBuffer + matSubImXY;
			
			%add averaged frame to correlation buffer
			if intTempAvgCounter == intTemporalAveraging
				intTempAvgCounter = 0;
				intBufferCounter = intBufferCounter + 1;
				matPixRespBuffer(:,:,intBufferCounter) = single(matPixRespAverageBuffer/intTemporalAveraging);
				matPixRespAverageBuffer(:) = 0;
			end
			
			if intBufferCounter == intBufferT
				fprintf('Processing %s%s... Now at %d%% [%s]; starting pixel-correlation analysis\n', sRec.strSession, sRec.sProcLib.strRecording, fracNow, getTime);
				
				intBufferCounter = 0;
				intCorrsCounter = intCorrsCounter + 1;
				% initialize matrix for current run
				vecSurrIdx = -intSpatialWindow:intSpatialWindow;
				intSpatWinSize = intSpatialWindow*2+1;
				intDivFac = (intBufferT-1);
				
				%% loop through center pixels
				matTempPixelCorrMat = zeros(intSubY,intSubX,intSpatWinSize,intSpatWinSize,'single');
				for intPixY = 1:intSubY
					for intPixX=1:intSubX
						vecCenterTrace = squeeze(matPixRespBuffer(intPixY,intPixX,:));
						
						%% loop through surround
						for intSurrIdxY = 1:intSpatWinSize
							intOffsetY = vecSurrIdx(intSurrIdxY);
							intSurrPixY = intPixY + intOffsetY;
							if intSurrPixY < 1 || intSurrPixY > intSubY,continue;end
							
							for intSurrIdxX = 1:intSpatWinSize
								intOffsetX = vecSurrIdx(intSurrIdxX);
								intSurrPixX = intPixX + intOffsetX;
								if intSurrPixX < 1 || intSurrPixX > intSubX,continue;end
								
								%check if already defined
								if matTempPixelCorrMat(intPixY,intPixX,intSurrIdxY,intSurrIdxX) ~= 0,continue;end
								
								%otherwise, calculate
								vecSurroundTrace = squeeze(matPixRespBuffer(intSurrPixY,intSurrPixX,:));
								sglCorr = sum(xzscore(vecCenterTrace,1).*xzscore(vecSurroundTrace,1))/intDivFac;
								matTempPixelCorrMat(intPixY,intPixX,intSurrIdxY,intSurrIdxX) = sglCorr;
								
								%calculate other location
								intNextY = intPixY + intOffsetY;
								intNextX = intPixX + intOffsetX;
								intNextSurrIdxY = intSpatialWindow-intOffsetY+1;
								intNextSurrIdxX = intSpatialWindow-intOffsetX+1;
								
								%forward assignment
								matTempPixelCorrMat(intNextY,intNextX,intNextSurrIdxY,intNextSurrIdxX) = sglCorr;
							end
						end
					end
				end
				
				%add values from current run to overall matrix
				matPixelCorrMat = matPixelCorrMat + matTempPixelCorrMat;
				
				fprintf('Finished pixel-correlation analysis [%s]\n', getTime);
			end
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
	%msg
	fprintf('Finished data collection... Calulating pixel responsiveness\n');
	
	%remove missing repetitions
	intRemove = find(squeeze(any(isnan(matPixelNormResp(1,1,:,:)),3)),1,'first');
	if ~isempty(intRemove)
		matPixelNormResp(:,:,:,intRemove:end) = [];
		matPixelStimResp(:,:,:,intRemove:end) = [];
	end
	intRemoveBase = find(squeeze(any(isnan(matPixelBaseResp(1,1,:,:)),3)),1,'first');
	if ~isempty(intRemoveBase)
		matPixelBaseResp(:,:,:,intRemoveBase:end) = [];
	end
	%z-score responses
	matPixelMeanNormResp = mean(matPixelNormResp,4);
	matPixelZScoreStimResp = nan(sRec.sProcLib.y,sRec.sProcLib.x,intStimTypes);
	
	[matBaseZ,matBaseMu,matBaseSigma] = zscore(squeeze(matPixelBaseResp),[],3); %#ok<ASGLU>
	for intStim=1:intStimTypes
		%get stim resp
		matThisResp = mean(squeeze(matPixelStimResp(:,:,intStim,:)),3);
		
		%z-score normalized to baseline
		matRespZ = (matThisResp-matBaseMu)./matBaseSigma;
		
		%put into matrix
		matPixelZScoreStimResp(:,:,intStim) = matRespZ;
	end
	
	%calculate overall pixel correlations
	matPixelCorrMat = matPixelCorrMat / intCorrsCounter;
	
	%calculate maximum visual responsiveness
	%max stim response in sd's above baseline
	matPixelMaxResponsiveness = max(matPixelZScoreStimResp,[],3);
	
	%calculate orientation selectivity
	%[(max mean normalized stim resp) - (min mean normalized stim resp)]/ [(max mean normalized stim resp) + (min mean normalized stim resp)]
	matMaxResp = max(matPixelMeanNormResp,[],3);
	matMinResp = min(matPixelMeanNormResp,[],3);
	matPixelSelectivity = (matMaxResp - matMinResp);% .* mean(cat(3,abs(matMaxResp),abs(matMinResp)),3);
	[indexList,indexLow,indexHigh] = getOutliers(matPixelSelectivity,10);
	matPixelSelectivity(indexList) = 0;
	
	%smooth
	ptrSmoothFilter = fspecial('disk', 2);
	matPixelMaxResponsiveness = imfilter(matPixelMaxResponsiveness,ptrSmoothFilter,'replicate');
	matPixelSelectivity = imfilter(matPixelSelectivity,ptrSmoothFilter,'replicate');
	
	%% save data into sRec
	sPixResp = struct;
	%sPixResp.matPixelStimResp = matPixelStimResp;
	%sPixResp.matPixelBaseResp = matPixelBaseResp;
	%sPixResp.matPixelNormResp = matPixelNormResp;
	sPixResp.strInfo1 = 'matPixelSelectivity is orientation selectivity defined as: [(max z-scored normalized stim resp) - (min z-scored normalized stim resp)]/ [(max z-scored normalized stim resp) + (min z-scored normalized stim resp)]';
	sPixResp.matPixelSelectivity = matPixelSelectivity;
	sPixResp.strInfo2 = 'matPixelMaxResponsiveness is visual responsivess defined as: number of sd''s that mean response is above baseline response (in baseline sd''s) [for stim type with highest value]';
	sPixResp.matPixelMaxResponsiveness = matPixelMaxResponsiveness;
	sPixResp.matPixelCorrMat = matPixelCorrMat;
	sRec.sPixResp = sPixResp;
	
	%save file
	if ~exist('strNewDir','var')
		strNewDir = [sRec.sMD.strMasterDir sRec.sMD.strImgTarget sRec.strSession filesep sRec.sProcLib.strRecording filesep];
	end
	strFile = [strNewDir sRec.strSession sRec.sProcLib.strRecording '_preproPixCorr.mat'];
	save([strFile],'sRec') %save temp
	fprintf('Saved recording structure to %s\n',strFile)
	
	%% run automatic cell detection
	if boolCalcPixCorr
		%get HQ image & subsample
		imProcOrig = sRec.imAverage.Overlay.ProcHQ;
		imProc = (imProcOrig(1:2:end,:,:) + imProcOrig(2:2:end,:,:))/2;
		imProc = (imProc(:,1:2:end,:) + imProc(:,2:2:end,:))/2;
		
		%run automatic detection
		[matAutoCellDetectMasks,matBelongsToSameAsNeighbor] = doCellDetection(matPixelCorrMat,imProc);
		
		% transform to original size & save data
		matBelongsToSameAsNeighbor = imresize(imnorm(matBelongsToSameAsNeighbor),2);
		matAutoCellDetectMasks = imresize(matAutoCellDetectMasks,2);
		
		%save to sRec
		sRec.sPixResp.matBelongsToSameAsNeighbor = matBelongsToSameAsNeighbor;
		sRec.sPixResp.matAutoCellDetectMasks = matAutoCellDetectMasks;
		
		%% save file
		strFile = [strNewDir sRec.strSession sRec.sProcLib.strRecording '_preproAutoDetected.mat'];
		save([strFile],'sRec') %save temp
		%movefile(strFile,[strFile '.backup'],'f');
		%movefile([strFile '1'],strFile,'f');
		fprintf('Saved automatic cell detection structure to %s\n',strFile)
	end
end


function [vecFramesAP, vecNumberAP, vecSpikes, vecExpFit, vecSpikeTimes] = doDetectSpikes(vec_dFoF,dblSamplingFreq,dblSpikeTau,intBlockSize,dblThresholdFactor)
	%doDetectSpikes Detects activation events in Ca2+ dF/F0 traces
	%   [vecFramesAP, vecNumberAP, vecSpikes, vecExpFit, vecSpikeTimes] = doDetectSpikes(vec_dFoF,dblSamplingFreq,dblSpikeTau,intBlockSize,dblThresholdFactor)
	%
	%	Version history:
	%	3.0 - July 20 2015
	%	Created by Jorrit Montijn
	%	4.0 - April 7 2020
	%	Added dynamic threshold and spike time vector output [by JM]
	
	%parameters & pre-allocate
	if ~exist('dblThresholdFactor','var') || isempty(dblThresholdFactor)
		dblThresholdFactor = 1;
	end
	
	if ~exist('intBlockSize','var') || isempty(intBlockSize)
		vecPrimes = primes(5000);
		intBlockSize = vecPrimes(end);
	elseif intBlockSize < 1 || round(intBlockSize) ~= intBlockSize
		error([mfilename 'E:SyntaxError'],'Syntax error; block size is incorrect');
	end
	
	%separate into blocks to avoid block edge intersections
	intOffset1 = 1;
	vecPrimes1000 = primes(1000);
	intOffset2 = vecPrimes1000(end);
	vecPrimes2000 = primes(2000);
	intOffset3 = vecPrimes2000(end);
	
	%get blocks with AE detections
	[vecSpikes1, vecExpFit1] = doDetectSpikeBlock(vec_dFoF(intOffset1:end),dblSamplingFreq,dblSpikeTau,intBlockSize,dblThresholdFactor);
	
	[vecSpikes2, vecExpFit2] = doDetectSpikeBlock(vec_dFoF(intOffset2:end),dblSamplingFreq,dblSpikeTau,intBlockSize,dblThresholdFactor);
	vecSpikes2 = [zeros(1,intOffset2-1) vecSpikes2];
	vecExpFit2 = [zeros(1,intOffset2-1) vecExpFit2];
	
	[vecSpikes3, vecExpFit3] = doDetectSpikeBlock(vec_dFoF(intOffset3:end),dblSamplingFreq,dblSpikeTau,intBlockSize,dblThresholdFactor);
	vecSpikes3 = [zeros(1,intOffset3-1) vecSpikes3];
	vecExpFit3 = [zeros(1,intOffset3-1) vecExpFit3];
	
	
	%combine vecSpikes
	vecSpikes = max(cat(1,vecSpikes1,vecSpikes2,vecSpikes3),[],1);
	
	%combine expfit
	vecExpFit = max(cat(1,vecExpFit1,vecExpFit2,vecExpFit3),[],1);
	
	%get frames
	vecFramesAP = find(vecSpikes>0);
	vecNumberAP = vecSpikes(vecSpikes>0);
	
	%% build spike times
	if nargout > 4
		%pre-allocate
		dblSpikeWindow = (1/dblSamplingFreq);
		vecUniqueSpikeNum = unique(vecNumberAP);
		vecSpikeTimes = nan(1,sum(vecNumberAP));
		intCounter = 1;
		%run through # of spikes per frame
		for intSpikesPerFrame=vecUniqueSpikeNum(:)'
			%get spikes and sort
			vecTheseSpikes= find(vecSpikes==intSpikesPerFrame);
			[vecSorted,vecOrderIdx] = sort(vecExpFit(vecTheseSpikes),'descend');
			
			%assign lag based on exp fit; higher is earlier
			intNumEvents = numel(vecOrderIdx);
			vecLag = ((1:intNumEvents)/intNumEvents)*dblSpikeWindow;
			%add original frame time
			vecCenterT = (((vecTheseSpikes(vecOrderIdx)-1)/dblSamplingFreq) + vecLag);
			%disperse spikes randomly in following window
			matRand = dblSpikeWindow*rand(intSpikesPerFrame,intNumEvents);
			%add dispersion to center time
			matSpikeTimes = bsxfun(@plus,matRand,vecCenterT);
			
			%assign spike times
			vecSpikeTimes(intCounter:(intCounter+(intNumEvents*intSpikesPerFrame)-1)) = matSpikeTimes(:);
			intCounter = intCounter + intNumEvents*intSpikesPerFrame;
		end
		
		%sort spike times
		vecSpikeTimes = sort(vecSpikeTimes,'ascend');
	end
end
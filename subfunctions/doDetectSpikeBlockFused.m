function [vecSpikes, vecExpFit] = doDetectSpikeBlockFused(vec_dFoF,dblSamplingFreq,dblSpikeTau,intBlockSize, dblThresholdFactor)
	%doDetectSpikeBlock Detects activation events in Ca2+ dF/F0 traces
	%   [vecFramesAP, vecNumberAP, vecSpikes, vecExpFit] = doDetectSpikeBlock(dFoF,dblSamplingFreq,dblSpikeTau,intBlockSize, dblThresholdFactor)
	%
	%	Version history:
	%	3.0 - January 29 2016
	%	Created by Jorrit Montijn
	%	4.0 - April 7 2020
	%	Added dynamic threshold
	%	4.1 - September 1 2020
	%	Reduced computation time and fixed some bugs
	
	%threshold factor
	if ~exist('dblThresholdFactor','var') || isempty(dblThresholdFactor)
		dblThresholdFactor = 1;
	end
	
	%parameters & pre-allocate
	if nargin < 4 || isempty(intBlockSize)
		intBlockSize = 997;
	end
	vecSpikes = nan(size(vec_dFoF));
	vecExpFit = nan(size(vec_dFoF));
	
	%get number of blocks (split by [intBlockSize] frames)
	intTotDur = length(vec_dFoF);
	intBlocks = ceil(intTotDur/intBlockSize);
	
	
	%split into blocks
	for intBlock=1:intBlocks
		%get data
		vecFrames = ((intBlock-1)*intBlockSize+1):min([(intBlock*intBlockSize) intTotDur]);
		vec_dFoF_block = vec_dFoF(vecFrames);
		intLength = numel(vecFrames);
		
		%append end to dF/F if last block
		if intBlock==intBlocks
			intAddFrames = min([100,numel(vec_dFoF_block)]);
			vec_dFoF_block((end+1):(end+intAddFrames)) = vec_dFoF_block(1:intAddFrames);
		end
		
		%perform spike detection on block
		
		% dFoF-criteria based transient detection to get initial guess
		vecTransients = getTransientGuess(vec_dFoF_block, dblSamplingFreq, dblSpikeTau, dblThresholdFactor);
		
		% loop for every block of transients to remove spikes with insufficient amplitude
		[apFrames, apSpikes, expFitB] = ...
			find_action_potentials_in_transients(vecTransients,vec_dFoF_block, dblSamplingFreq, dblSpikeTau, dblThresholdFactor);
		
		%vectorized recombination
		vecSpikesB = zeros(size(vec_dFoF_block));
		vecSpikesB(apFrames) = apSpikes;
		
		%assign to outpit
		vecSpikes(vecFrames) = vecSpikesB(1:intLength);
		vecExpFit(vecFrames) = expFitB(1:intLength);
	end
end

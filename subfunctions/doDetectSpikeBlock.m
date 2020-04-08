function [vecSpikes, vecExpFit] = doDetectSpikeBlock(vec_dFoF,dblSamplingFreq,dblSpikeTau,intBlockSize, dblThresholdFactor)
	%doDetectSpikeBlock Detects activation events in Ca2+ dF/F0 traces
	%   [vecFramesAP, vecNumberAP, vecSpikes, vecExpFit] = doDetectSpikeBlock(dFoF,dblSamplingFreq,dblSpikeTau,intBlockSize, dblThresholdFactor)
	%
	%	Version history:
	%	3.0 - January 29 2016
	%	Created by Jorrit Montijn
	%	4.0 - April 7 2020
	%	Added dynamic threshold
	
	%threshold factor
	if ~exist('dblThresholdFactor','var') || isempty(dblThresholdFactor)
		dblThresholdFactor = 1;
	end
	
	%parameters & pre-allocate
	if nargin < 4 || isempty(intBlockSize)
		intBlockSize = 2500;
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
			intAddFrames = 100;
			vec_dFoF_block((end+1):(end+intAddFrames)) = vec_dFoF_block(1:intAddFrames);
		end
		
		%perform spike detection on block
		
		% dFoF-criteria based transient detection to get initial guess
		vecTransients = getTransientGuess(vec_dFoF_block, dblSamplingFreq, dblSpikeTau, dblThresholdFactor);
		
		% separate into blocks with connected transients
		[sep_transients, sep_dFoF, sep_start, sep_stop] = ...
			separate_transients( vecTransients, vec_dFoF_block, dblSamplingFreq, dblSpikeTau );
		
		% loop for every block of transients to remove spikes with insufficient amplitude
		sep_apFrames = cell(1,length(sep_transients)); 
		sep_apSpikes = cell(1,length(sep_transients)); 
		sep_expFit = cell(1,length(sep_transients));
		for b = 1:length(sep_transients);
			[sep_apFrames{b}, sep_apSpikes{b}, sep_expFit{b}] = ...
				find_action_potentials_in_transients( sep_transients{b}, sep_dFoF{b}, dblSamplingFreq, dblSpikeTau, dblThresholdFactor);
		end
		
		% recombine into one trace and one list of AP's
		[vecSpikesB, expFitB] = doRecombineTransients(vec_dFoF_block, sep_apFrames, sep_apSpikes, sep_expFit, sep_start, sep_stop);
		
		%assign to outpit
		vecSpikes(vecFrames) = vecSpikesB(1:intLength);
		vecExpFit(vecFrames) = expFitB(1:intLength);
	end
end

function [sep_transients, sep_dFoF, sep_start, sep_stop] = ...
		separate_transients( vecTransients, vec_dFoF_block, dblSamplingFreq, dblSpikeTau )
	%dt
	dt = 1/dblSamplingFreq ;
	
	% number of transients
	m = length(vecTransients) ;
	intLength = 50;
	
	% put all transients in one array
	ExpMat = zeros(m, length(vec_dFoF_block)) ;
	for j = 1:m
		ExpMat( j, vecTransients(j):vecTransients(j)+intLength ) = ...
			exp( (-(0:1:intLength) .* dt) ./ dblSpikeTau ) ;
	end
	transarray = sum(ExpMat);
	
	% extract blocks of connected non-zero timepoints
	blnr = 0;
	blstart = 1;
	blstop = 0;
	sep_transients = cell(0);
	sep_dFoF = cell(0);
	intBlockMax = 1000;
	sep_start = nan(1,intBlockMax);
	sep_stop = nan(1,intBlockMax);
	
	
	for t = 1:length(transarray);
		if transarray(t) == 0
			% part of data that doesn't belong to a block is detected
			
			if blstart < blstop
				% end of block detected, add transients and dFoF to output
				% variables
				blnr = blnr + 1;
				
				%pre-allocate in chunks of 1000
				if blnr > intBlockMax
					intBlockMax = intBlockMax + 1000;
					sep_start = [sep_start nan(1,1000)];
					sep_stop = [sep_stop nan(1,1000)];
				end
				
				% cut out trace of dFoF
				sep_dFoF{blnr} = vec_dFoF_block(blstart:blstop);
				
				% find transients that go in trace
				trns = vecTransients(vecTransients >= blstart & vecTransients < blstop);
				
				% correct for changed starting frame
				trns = trns - (blstart-1);
				sep_transients{blnr} = trns;
				
				% update block start and stop
				sep_start(blnr) = blstart;
				sep_stop(blnr) = blstop;
				blstart = 1;
				blstop = 0;
				
			end
		else
			% part of data belongs to block
			
			if blstart > blstop
				% start of block detected
				blstart = t;
				blstop = t;
			else
				% next data point in block
				blstop = t;
			end
		end
	end
	sep_start = sep_start(~isnan(sep_start));
	sep_stop = sep_stop(~isnan(sep_stop));
end

function [apFrames, apSpikes, expFit] = ...
		find_action_potentials_in_transients( vecSepTransients, vecSep_dFoF, dblSamplingFreq, dblSpikeTau, dblThresholdFactor)
	
	intLength = 50;
	
	if size(vecSep_dFoF,1) == 1
		vecSep_dFoF = vecSep_dFoF';
	end
	dt = 1/dblSamplingFreq ;
	
	% number of transients
	m = length(vecSepTransients) ;
	
	% matrix with one row per transient and one exponential on each row,
	% starting at the location of the transient
	ExpMat = zeros(m, length(vecSep_dFoF)) ;
	for j = 1:m
		ExpMat( j, vecSepTransients(j):vecSepTransients(j)+intLength ) = ...
			exp( (-(0:1:intLength) .* dt) ./ dblSpikeTau ) ;
	end
	
	% find weight of transient at each row that gives best fit to the dFoF
	trBelow01 = 1;
	while trBelow01
		
		% fit transients with matrix of exponentials
		%         h = lsqr(ExpMat', dFoF, 0.001, 50);
		h = ExpMat'\vecSep_dFoF;
		
		% find exponential with the smallest height
		minh = find( h == min(h) );
		if length(minh) > 1
			minh = minh(1);
		end
		
		% check if minh < 0.1
		if h(minh) < (0.1*dblThresholdFactor)
			% remove transient and continue loop
			ExpMat(minh,:) = [];
			vecSepTransients(minh) = [];
		else
			% no transients are smaller that 0.1, so we're done
			trBelow01 = 0;
		end
		
	end
	
	% set output variables
	apFrames = vecSepTransients;
	apSpikes = round(((h*100)/10)/dblThresholdFactor);
	expFit = ExpMat'*h;
	
end
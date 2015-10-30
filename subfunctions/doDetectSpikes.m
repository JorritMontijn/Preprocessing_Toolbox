function [apFrames, apSpikes, vecSpikes, expFit] = doDetectSpikes(dFoF,dblSamplingFreq,dblTau)
	%doDetectSpikes Detects activation events in Ca2+ dF/F0 traces
	%   [apFrames, apSpikes, vecSpikes, expFit] = doDetectSpikes(dFoF,dblSamplingFreq,dblTau)
	%
	%	Version history:
	%	2.0 - July 20 2015
	%	Created by Jorrit Montijn
	
	%parameters & pre-allocate
	intBlockSize = 20000;
	apFrames = [];
	apSpikes = [];
	vecSpikes = nan(size(dFoF));
	expFit = nan(size(dFoF));
	
	%get number of blocks (split by [intBlockSize] frames)
	intTotDur = length(dFoF);
	intBlocks = ceil(intTotDur/20000);
	
	
	%split into blocks
	for intBlock=1:intBlocks
		%get data
		vecFrames = ((intBlock-1)*intBlockSize+1):min([(intBlock*intBlockSize) intTotDur]);
		vec_dFoF_block = dFoF(vecFrames);
		
		%perform spike detection on block
		
		% dFoF-criteria based transient detection to get initial guess
		[transients] = getTransientGuess(vec_dFoF_block, dblSamplingFreq, dblTau);
		
		% separate into blocks with connected transients
		[sep_transients, sep_dFoF, sep_start, sep_stop] = ...
			separate_transients( transients, vec_dFoF_block, dblSamplingFreq, dblTau );
		
		% loop for every block of transients to remove spikes with insufficient amplitude
		sep_apFrames = cell(0); sep_apSpikes = cell(0); sep_expFit = cell(0);
		for b = 1:length(sep_transients);
			[sep_apFrames{b}, sep_apSpikes{b}, sep_expFit{b}] = ...
				find_action_potentials_in_transients( sep_transients{b}, sep_dFoF{b}, dblSamplingFreq, dblTau );
		end
		
		% recombine into one trace and one list of AP's
		[apFramesB, apSpikesB, vecSpikesB, expFitB] = doRecombineTransients(vec_dFoF_block, sep_apFrames, sep_apSpikes, sep_expFit, sep_start, sep_stop);
		
		
		%assign to outpit
		apFrames = [apFrames apFramesB]; %#ok<AGROW>
		apSpikes = [apSpikes apSpikesB]; %#ok<AGROW>
		vecSpikes(vecFrames) = vecSpikesB;
		expFit(vecFrames) = expFitB;
	end
end

function [sep_transients, sep_dFoF, sep_start, sep_stop] = ...
		separate_transients( transients, dFoF, samplingFreq, tau )
	dt = 1/samplingFreq ;
	
	% number of transients
	m = length(transients) ;
	intLength = 50;
	
	% put all transients in one array
	ExpMat = zeros(m, length(dFoF)) ;
	for j = 1:m
		ExpMat( j, transients(j):transients(j)+intLength ) = ...
			exp( (-(0:1:intLength) .* dt) ./ tau ) ;
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
				sep_dFoF{blnr} = dFoF(blstart:blstop);
				
				% find transients that go in trace
				trns = transients(transients >= blstart & transients < blstop);
				
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
		find_action_potentials_in_transients( transients, dFoF, samplingFreq, tau )
	
	intLength = 50;
	
	if size(dFoF,1) == 1
		dFoF = dFoF';
	end
	dt = 1/samplingFreq ;
	
	% number of transients
	m = length(transients) ;
	
	% matrix with one row per transient and one exponential on each row,
	% starting at the location of the transient
	ExpMat = zeros(m, length(dFoF)) ;
	for j = 1:m
		ExpMat( j, transients(j):transients(j)+intLength ) = ...
			exp( (-(0:1:intLength) .* dt) ./ tau ) ;
	end
	
	% find weight of transient at each row that gives best fit to the dFoF
	trBelow01 = 1;
	while trBelow01
		
		% fit transients with matrix of exponentials
		%         h = lsqr(ExpMat', dFoF, 0.001, 50);
		h = ExpMat'\dFoF;
		
		% find exponential with the smallest height
		minh = find( h == min(h) );
		if length(minh) > 1
			minh = minh(1);
		end
		
		% check if minh < 0.1
		if h(minh) < 0.1
			% remove transient and continue loop
			ExpMat(minh,:) = [];
			transients(minh) = [];
		else
			% no transients are smaller that 0.1, so we're done
			trBelow01 = 0;
		end
		
	end
	
	% set output variables
	apFrames = transients;
	apSpikes = round((h*100)/9.75);
	expFit = ExpMat'*h;
	
end
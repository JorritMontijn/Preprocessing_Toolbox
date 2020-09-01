
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


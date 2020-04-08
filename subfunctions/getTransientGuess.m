function vecTransients = getTransientGuess(vec_dFoF, dblSamplingFreq, dblSpikeTau, dblThresholdFactor)
	%getTransientGuess Fast heuristic search for potential transients (based on Greenberg et al.)
	%   Syntax: vecTransients = getTransientGuess(vec_dFoF, dblSamplingFreq, dblSpikeTau, dblThresholdFactor)
	
	if ~exist('dblThresholdFactor','var') || isempty(dblThresholdFactor)
		dblThresholdFactor = 1;
	end
	
	if size(vec_dFoF,1) == 1
		vec_dFoF = vec_dFoF';
	end
	intBlockMax = 1000;
	vecTransients = nan(1,intBlockMax);
	
	t = 0 ;
	dt = 1/dblSamplingFreq ;
	Y = exp( (-(0:1:17) .* dt) ./ dblSpikeTau ) ;
	
	%run vectorized heuristic
	vecSelect = 3:(length(vec_dFoF)-17);
	indA = vec_dFoF(vecSelect) > (0.06*dblThresholdFactor) ;
	indB = (vec_dFoF(vecSelect) - vec_dFoF(vecSelect-1)) > (0.02*dblThresholdFactor);
	indC = (vec_dFoF(vecSelect) - vec_dFoF(vecSelect-2)) > (0.008*dblThresholdFactor);
	indD = (vec_dFoF(vecSelect+1) - vec_dFoF(vecSelect-1)) > (-0.03/dblThresholdFactor);
	
	for i=find(indA & indB & indC & indD)'
		eqE = (vec_dFoF(i:i+17)' * Y') / sqrt( sum(Y.^2) ) ;
		if eqE > (0.08*dblThresholdFactor)
			t = t + 1 ;
			if t > intBlockMax
				intBlockMax = intBlockMax + 1000;
				vecTransients = [vecTransients nan(1,1000)];
			end
			vecTransients(t) = i ;
		end
	end
	
	vecTransients = vecTransients(~isnan(vecTransients));
end


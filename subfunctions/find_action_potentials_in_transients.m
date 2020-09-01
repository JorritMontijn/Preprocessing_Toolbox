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
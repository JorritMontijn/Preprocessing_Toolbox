function [transients] = getTransientGuess(dFoF, samplingFreq, tau)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	
	if size(dFoF,1) == 1
		dFoF = dFoF';
	end
	intBlockMax = 1000;
	transients = nan(1,intBlockMax);
	
	t = 0 ;
	dt = 1/samplingFreq ;
	Y = exp( (-(0:1:17) .* dt) ./ tau ) ;
	
	
	intType =2;
	if intType == 1
		% check which frames meet the initial criteria
		
		for i = 3:(length(dFoF)-17)
			eqA = dFoF(i) ;
			eqB = dFoF(i) - dFoF(i-1) ;
			eqC = dFoF(i) - dFoF(i-2) ;
			eqD = dFoF(i+1) - dFoF(i-1) ;
			eqE = (dFoF(i:i+17)' * Y') / sqrt( sum(Y.^2) ) ;
			
			if eqA > 0.06 && eqB > 0.02 && eqC > 0.008 && eqD > -0.03 && eqE > 0.08
				t = t + 1 ;
				if t > intBlockMax
					intBlockMax = intBlockMax + 1000;
					transients = [transients nan(1,1000)];
				end
				transients(t) = i ;
			end
		end
		
		
	else
		vecSelect = 3:(length(dFoF)-17);
		indA = dFoF(vecSelect) > 0.06 ;
		indB = (dFoF(vecSelect) - dFoF(vecSelect-1)) > 0.02;
		indC = (dFoF(vecSelect) - dFoF(vecSelect-2)) > 0.008;
		indD = (dFoF(vecSelect+1) - dFoF(vecSelect-1)) > -0.03;
		
		for i=find(indA & indB & indC & indD)'
			eqE = (dFoF(i:i+17)' * Y') / sqrt( sum(Y.^2) ) ;
			if eqE > 0.08
				t = t + 1 ;
				if t > intBlockMax
					intBlockMax = intBlockMax + 1000;
					transients = [transients nan(1,1000)];
				end
				transients(t) = i ;
			end
		end
		
	end
	transients = transients(~isnan(transients));
end


function [vecSpikes, expFit] = doRecombineTransients(dFoF, sep_apFrames, sep_apSpikes, sep_expFit, sep_start, sep_stop)
	%doRecombineTransients Recombines transients
	%   Syntax: [vecSpikes, expFit] = doRecombineTransients(dFoF, sep_apFrames, sep_apSpikes, sep_expFit, sep_start, sep_stop)
	
	%vectorized recombination
	expFit = zeros(size(dFoF));
	vecSpikes = zeros(size(dFoF));
	for b = 1:length(sep_apFrames)
		vecSpikes(sep_apFrames{b} + (sep_start(b)-1)) = sep_apSpikes{b};
		expFit(sep_start(b):sep_stop(b)) = sep_expFit{b};
	end
end


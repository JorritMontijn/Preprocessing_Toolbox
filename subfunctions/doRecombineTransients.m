function [apFrames, apSpikes, vecSpikes, expFit] = doRecombineTransients(dFoF, sep_apFrames, sep_apSpikes, sep_expFit, sep_start, sep_stop)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	intType = 2;
	if intType == 1
		apFrames = [];
		apSpikes = [];
		expFit = zeros(size(dFoF));
		vecSpikes = zeros(size(dFoF));
		
		% recombine multiple lists of ap's
		for b = 1:length(sep_apFrames)
			
			% recombine apFrames
			apFrames(end+1:end+length(sep_apFrames{b})) = sep_apFrames{b} + (sep_start(b)-1);
			
			% recombine apSpikes
			apSpikes(end+1:end+length(sep_apSpikes{b})) = sep_apSpikes{b};
			
			% recombine expFit
			expFit( sep_start(b):sep_stop(b),1 ) = sep_expFit{b};
			
		end
		vecSpikes(apFrames) = apSpikes;
	else
		expFit = zeros(size(dFoF));
		vecSpikes = zeros(size(dFoF));
		for b = 1:length(sep_apFrames)
			vecSpikes(sep_apFrames{b} + (sep_start(b)-1)) = sep_apSpikes{b};
			expFit(sep_start(b):sep_stop(b)) = sep_expFit{b};
		end
		indSpikes = vecSpikes>0;
		apSpikes = vecSpikes(indSpikes);
		apFrames = find(indSpikes);
	end
end


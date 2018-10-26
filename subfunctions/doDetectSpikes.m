function [apFrames, apSpikes, vecSpikes, expFit] = doDetectSpikes(dFoF,dblSamplingFreq,dblTau,intBlockSize)
	%doDetectSpikes Detects activation events in Ca2+ dF/F0 traces
	%   [apFrames, apSpikes, vecSpikes, expFit] = doDetectSpikes(dFoF,dblSamplingFreq,dblTau)
	%
	%	Version history:
	%	3.0 - July 20 2015
	%	Created by Jorrit Montijn
	
	%parameters & pre-allocate
	if nargin < 4 || isempty(intBlockSize)
		intBlockSize = 1000;
	end
	
	if intBlockSize < 5000
		intOffset1 = 1;
		intOffset2 = 200;
		intOffset3 = 400;
		
		%get blocks with AE detections
		[apFrames1, apSpikes1, vecSpikes1, expFit1] = doDetectSpikeBlock(dFoF(intOffset1:end),dblSamplingFreq,dblTau,intBlockSize);
		
		[apFrames2, apSpikes2, vecSpikes2, expFit2] = doDetectSpikeBlock(dFoF(intOffset2:end),dblSamplingFreq,dblTau,intBlockSize);
		vecSpikes2 = [zeros(1,intOffset2-1) vecSpikes2];
		expFit2 = [zeros(1,intOffset2-1) expFit2];
		
		[apFrames3, apSpikes3, vecSpikes3, expFit3] = doDetectSpikeBlock(dFoF(intOffset3:end),dblSamplingFreq,dblTau,intBlockSize);
		vecSpikes3 = [zeros(1,intOffset3-1) vecSpikes3];
		expFit3 = [zeros(1,intOffset3-1) expFit3];
		
		
		%combine vecSpikes
		vecSpikes = (vecSpikes1 + vecSpikes2 + vecSpikes3) / 3;
		vecSpikes(1:200) = vecSpikes1(1:200);
		vecSpikes(201:400) = (vecSpikes2(201:400) + vecSpikes3(201:400))/2;
		vecSpikes = ceil(vecSpikes);
		
		%combine expfit
		expFit = (expFit1 + expFit2 + expFit3) / 3;
		expFit(1:200) = expFit1(1:200);
		expFit(201:400) = (expFit1(201:400) + expFit2(201:400))/2;
		
		%get frames
		apFrames = find(vecSpikes>0);
		apSpikes = vecSpikes(vecSpikes>0);
		
	else
		%get AE detections
		[apFrames, apSpikes, vecSpikes, expFit] = doDetectSpikeBlock(dFoF,dblSamplingFreq,dblTau,intBlockSize);
	end
end
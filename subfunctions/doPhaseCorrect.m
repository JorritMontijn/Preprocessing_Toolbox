function imCorrected = doPhaseCorrect(imCorrected,intShift)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	imShift = circshift(imCorrected,[0 intShift]);
	imCorrected(2:2:end,:) = imShift(2:2:end,:);
end


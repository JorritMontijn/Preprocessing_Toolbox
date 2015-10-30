function [dblPhaseCorrect,output] = doPhaseRegistration(imUncorr)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	matUncorr = im2double(imUncorr);
	if ndims(matUncorr) == 3
		matUncorr = sum(matUncorr,3);
	end
	
	imOdd = matUncorr(1:2:end,:);
	imEven = matUncorr(2:2:end,:);
	
	output = dftregistration(fft2(imOdd),fft2(imEven),1);
	dblPhaseCorrect = output(4);
end


function matIm = doRemSaturated(matIm)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	if max(matIm(:)) == intmax('uint8')
		matSatIndex = matIm == intmax('uint8');
		matIm = doReplaceBySurround(matIm,matSatIndex);
	end
end
function matOut = doReplaceBySurround(matIm,matReplaceIndex)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	[maxY,maxX] = size(matIm);
	[vecY,vecX] = find(matReplaceIndex);
	matIm(matReplaceIndex) = nan;
	matOut = matIm;
	for intPix=1:length(vecX)
		intX=vecX(intPix);
		intY=vecY(intPix);
		
		%define masks
		maskX=max(min([intX-1 intX intX+1 ...
			intX-1 intX+1 ...
			intX-1 intX intX+1],maxX),1);
		
		maskY=max(min([intY-1 intY-1 intY-1 ...
			intY intY ...
			intY+1 intY+1 intY+1],maxY),1);
		
		%get surrounding pixels
		vecSurround = diag(matIm(maskY,maskX));
		
		%get max surround and assign to location
		matOut(intY,intX) = nanmean(vecSurround);
	end
end


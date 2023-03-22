function matIm = doRemSaturated(matIm)
	%doRemSaturated Remove saturated pixels in 2D matrix
	%   matIm = doRemSaturated(matIm)
	
	if isfloat(matIm)
		dblMax = 1;
	else
		dblMax = intmax(class(matIm));
	end
	if max(matIm(:)) == dblMax
		matSatIndex = matIm == dblMax;
		matIm = doReplaceBySurround(matIm,matSatIndex);
	end
end
function matOut = doReplaceBySurround(matIm,matReplaceIndex)
	%doReplaceBySurround Replace entries in 2D matrix by max of surround
	%   matOut = doReplaceBySurround(matIm,matReplaceIndex)
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
		vecSurround = flat(matIm(maskY,maskX));
		
		%get max surround and assign to location
		matOut(intY,intX) = nanmean(vecSurround);
	end
end


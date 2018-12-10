function [matMeanCorrOverlap,matSDCorrOverlap] = DC_ACD_GetCrossCorrPix(matPixelCorrMat,matOverlapMask,matObjectMask)
	%get internal correlations
	
	%subsample object mask
	matSubMaskObject = matObjectMask(1:2:end,:) + matObjectMask(2:2:end,:);
	matSubMaskObject = matSubMaskObject(:,1:2:end) + matSubMaskObject(:,2:2:end);
	matSubMaskObject = matSubMaskObject > 1;
	
	%subsample overlap mask
	matSubMaskOverlap = matOverlapMask(1:2:end,:) + matOverlapMask(2:2:end,:);
	matSubMaskOverlap = matSubMaskOverlap(:,1:2:end) + matSubMaskOverlap(:,2:2:end);
	matSubMaskOverlap = matSubMaskOverlap > 1;
	matSubMeanCorrOverlap = zeros(size(matSubMaskObject));
	matSubSDCorrOverlap = zeros(size(matSubMaskObject));
	
	%get mask pixels
	[intSubSizeY,intSubSizeX] = size(matSubMaskOverlap);
	[vecRow,vecCol] = find(matSubMaskOverlap);
	intOverlapPixels = numel(vecRow);
	intSpatWin = ((size(matPixelCorrMat,3)-1)/2);
	vecSurround = -intSpatWin:intSpatWin;
	for intPix=1:intOverlapPixels
		intY = vecRow(intPix);
		intX = vecCol(intPix);
		matCorr = squeeze(matPixelCorrMat(intY,intX,:,:));
		
		%get subvectors
		vecY = vecSurround+intY;
		indRemY = vecY<1 | vecY>intSubSizeY;
		vecY(indRemY) = [];
		vecX = vecSurround+intX;
		indRemX = vecX<1 | vecX>intSubSizeX;
		vecX(indRemX) = [];
		
		%get submask
		matObjectMaskSurround = matSubMaskObject(vecY,vecX);
		matSubCorr = matCorr(~indRemY,~indRemX);
		vecCorrVals = matSubCorr(matObjectMaskSurround);
		
		%save
		matSubMeanCorrOverlap(intY,intX) = xmean(vecCorrVals,1);
		matSubSDCorrOverlap(intY,intX) = xstd(vecCorrVals,1);
	end
	
	%expand
	matMeanCorrOverlap = zeros(size(matOverlapMask));
	matMeanCorrOverlap(1:2:end,1:2:end) = matSubMeanCorrOverlap;
	matMeanCorrOverlap(1:2:end,2:2:end) = matSubMeanCorrOverlap;
	matMeanCorrOverlap(2:2:end,1:2:end) = matSubMeanCorrOverlap;
	matMeanCorrOverlap(2:2:end,2:2:end) = matSubMeanCorrOverlap;
	
	matSDCorrOverlap = zeros(size(matOverlapMask));
	matSDCorrOverlap(1:2:end,1:2:end) = matSubSDCorrOverlap;
	matSDCorrOverlap(1:2:end,2:2:end) = matSubSDCorrOverlap;
	matSDCorrOverlap(2:2:end,1:2:end) = matSubSDCorrOverlap;
	matSDCorrOverlap(2:2:end,2:2:end) = matSubSDCorrOverlap;
end
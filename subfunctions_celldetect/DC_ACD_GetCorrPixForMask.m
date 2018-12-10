function [dblMeanCorr,dblSDCorr,vecMeanCorr,vecSDCorr] = DC_ACD_GetCorrPixForMask(matMask,matPixelCorrMat)
	%get internal correlations
	
	%subsample mask
	matSubMask = matMask(1:2:end,:) + matMask(2:2:end,:);
	matSubMask = matSubMask(:,1:2:end) + matSubMask(:,2:2:end);
	matSubMask = matSubMask > 1;
	
	%get mask pixels
	[intSubSizeY,intSubSizeX] = size(matSubMask);
	[vecRow,vecCol] = find(matSubMask);
	intSubPixels = numel(vecRow);
	vecMeanCorr = nan(1,intSubPixels);
	vecSDCorr = nan(1,intSubPixels);
	intSpatWin = ((size(matPixelCorrMat,3)-1)/2);
	vecSurround = -intSpatWin:intSpatWin;
	for intPix=1:intSubPixels
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
		matSurrSubMask = matSubMask(vecY,vecX);
		matSubCorr = matCorr(~indRemY,~indRemX);
		vecCorrVals = matSubCorr(matSurrSubMask);
		
		%save
		vecMeanCorr(intPix) = xmean(vecCorrVals,1);
		vecSDCorr(intPix) = xstd(vecCorrVals,1);
	end
	
	%get scalar outputs
	dblMeanCorr = mean(vecMeanCorr);
	dblSDCorr = std(vecMeanCorr);
end
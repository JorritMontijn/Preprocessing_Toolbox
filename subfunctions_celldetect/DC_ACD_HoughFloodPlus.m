function matMask = DC_ACD_HoughFloodPlus(dblX, dblY, im1D)
	%% get images
	%global sFig;
	dblX =  9.3726;
	dblY = 39.5811;
	global sRec;
	if ~exist('im1D','var') || isempty(im1D)
		im1D = sRec.sPixResp.matFloodBorders;
	end
	[intMaxY,intMaxX] = size(im1D);
	
	%% get surround flood
	vecFilt = normpdf(-2:2,0,1);
	vecFilt = vecFilt ./ sum(vecFilt);
	
	intX = round(dblX);
	intY = round(dblY);
	vecPix = -50:50;
	vecX = vecPix+intX;
	vecX(vecX<1 | vecX>intMaxX)=[];
	vecY = vecPix+intY;
	vecY(vecY<1 | vecY>intMaxY)=[];
	matSubFlood = im1D(vecY,vecX);
	
	%find potential disk centers
	intSubCenterX = find(vecX==intX);
	intSubCenterY = find(vecY==intY);
	[matMergedCentroids,vecMergedRadii] = DC_ACD_GetDisks(matSubFlood,5);
	if isempty(matMergedCentroids)
		dblUseX = 0;
		dblUseY = 0;
		dblRadius = [];
	else
		vecDistsToCenter = sqrt((matMergedCentroids(:,1) - intSubCenterX).^2 + (matMergedCentroids(:,2) - intSubCenterY).^2);
		[dummy,intCentral] = min(vecDistsToCenter);
		vecCenter = matMergedCentroids(intCentral,:); %xy
		dblUseX = vecCenter(1)-intSubCenterX;
		dblUseY = vecCenter(2)-intSubCenterY;
		dblRadius = vecMergedRadii(intCentral);
	end
	
	%get distance to center & angle
	[matGridX,matGridY]=meshgrid(vecX-intX,vecY-intY);
	matDist = sqrt((matGridX-dblUseX).^2 + (matGridY-dblUseY).^2);
	matTheta = atan2(matGridX-dblUseX,matGridY-dblUseY);
	
	%get annulus
	if ~exist('dblRadius','var') || isempty(dblRadius) || dblRadius < 5
		[dummy,vecPeakVals] = makeBins(matDist(:),matSubFlood(:),1:1:(max(matDist)+1)); %#ok<*ASGLU>
		[dummy,vecRadii,vecWidths,vecHeights] = findpeaks(vecPeakVals);
		[dummy,intUsePeak] = max(vecHeights);
		dblRadius = vecRadii(intUsePeak)-2;
		if dblRadius < 5 || dblRadius > 30 %try filter
			vecPeakVals = imfilter(vecPeakVals,vecFilt,'replicate');
			[dummy,vecRadii,vecWidths,vecHeights] = findpeaks(vecPeakVals);
			[dummy,intUsePeak] = max(vecHeights);
			dblRadius = vecRadii(intUsePeak)-2;
		end
		dblWidth = max(vecWidths(intUsePeak),dblRadius/4);
		dblMinRadius = dblRadius-dblWidth/2;
		dblMaxRadius = dblRadius+dblWidth/2;
	else
		dblWidth = dblRadius;
		dblMinRadius = dblRadius-dblWidth/4;
		dblMaxRadius = dblRadius+dblWidth/2;
	end
	
	%% split into polar regions
	intRegions = 12;
	vecDodecantRadii = dblRadius*ones(1,intRegions);
	vecThetaSpaceBorders = linspace(-pi,pi,intRegions+1);
	dblThetaStep = unique(roundi(diff(vecThetaSpaceBorders),10));
	vecThetaSpaceMeans = vecThetaSpaceBorders(2:end)-dblThetaStep/2;
	[dummy,dummy,dummy,dummy,cellIDs] = makeBins(matTheta(:),matSubFlood(:),vecThetaSpaceBorders);
	
	for intRegion=1:numel(cellIDs)
		indPixels = cellIDs{intRegion};
		vecTheseDists = matDist(indPixels);
		vecTheseVals = matSubFlood(indPixels);
		[dummy,vecMeanOfRho] = makeBins(vecTheseDists,vecTheseVals,1:1:(max(vecTheseDists)+1));
		vecMeanOfRho(isnan(vecMeanOfRho)) = 0;
		vecMeanOfRho = conv(vecMeanOfRho,vecFilt,'same');
		[dummy,vecRadiiRho] = findpeaks(vecMeanOfRho);
		if isempty(vecRadiiRho) || vecRadiiRho(1) < dblMinRadius || vecRadiiRho(1) > dblMaxRadius
			vecDodecantRadii(intRegion) = dblRadius-1;
		else
			vecDodecantRadii(intRegion) = vecRadiiRho(1);
		end
	end
	
	%% transform polygons to mask
	[vecObjectY,vecObjectX] = pol2cart(vecThetaSpaceMeans,vecDodecantRadii);
	
	%get mask
	matMask = poly2mask(vecObjectX+intX+dblUseX, vecObjectY+intY+dblUseY, sRec.sProcLib.x, sRec.sProcLib.y );
end
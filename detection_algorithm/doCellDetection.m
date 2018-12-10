function [matAutoCellDetectMasks,matBelongsToSameAsNeighbor] = doCellDetection(matPixelCorrMat,imProc)
	%doCellDetection Performs automatic cell detection
	%   Syntax: [matAutoCellDetectMasks,matBelongsToSameAsNeighbor] = doCellDetection(matPixelCorrMat,imProc)
	%
	%   Input:
	%	- matPixelCorrMat, pre-processed pix-wise correlations [X x Y x SurrX x SurrY]
	%	- imProc, average image [X x Y]
	%
	%	Output:
	%	- matAutoCellDetectMasks: [X x Y x ROI] logical matrix with masks
	%	- matBelongsToSameAsNeighbor [X x Y] double matrix with neighbour-similarity
	%
	%Additional information:
	% 
	
	%% outline
	%pre-allocate cell array for all pixels with cell array list of
	%pixels belonging to same cluster
	%cellClust{intPixel}{intOtherPixel}(intClusterMember)
	
	%go through all pixels
	
	%build pixel-wise distance matrix from correlations
	
	%cluster pixels in this neighborhood
	
	%threshold and save cluster (if any) to {intOtherPixel} for all
	%pixels that share its cluster
	
	%% set switches
	boolPlot = false;
	boolFast = false;
	
	%% calculate surround window
	matMean = squeeze(xmean(xmean(matPixelCorrMat,1),2));
	intSpatWin = size(matMean,1);
	intSubSW = floor(intSpatWin/2);
	intSubX = size(matPixelCorrMat,2);
	intSubY = size(matPixelCorrMat,1);
	[vecSurrX,vecSurrY]=meshgrid((-intSubSW):(intSubSW),(-intSubSW):(intSubSW));
	vecSurrX = vecSurrX(:);
	vecSurrY = vecSurrY(:);
	intSurrPixTot = numel(vecSurrX);
	intThisPixIdx = ceil(intSurrPixTot/2);
	
	%% build pixel lookup table
	[matX,matY] = meshgrid(1:intSubX,1:intSubY);
	vecPixIdx = 1:numel(matX);
	vecPixX = matX(:);
	vecPixY = matY(:);
	intPixTot = numel(vecPixIdx);
	
	%% pre-allocate cell array
	cellClust = cell(1,intPixTot);
	
	%% build aggregate map to select potential regions of interest
	intLastX = 0;
	for intPixelIdx=1:intPixTot
		intY = vecPixY(intPixelIdx);
		if boolFast&&mod(intY,3) ~= 0,continue;end %skip y-pixels
		intX = vecPixX(intPixelIdx);
		if boolFast&&mod(intX,3) ~= 0,continue;end %skip x-pixels
		if intX ~= intLastX,fprintf('Processing x=%d [%s]\n',intX,getTime);intLastX=intX;end
		
		%build pixel-wise distance matrix from correlations
		matCorr = nan(intSpatWin^2,intSpatWin^2);
		
		%% go through surround
		for intSurrPix=1:numel(vecSurrX)
			intSurrX = vecSurrX(intSurrPix)+intX;
			if intSurrX < 1 || intSurrX > intSubX,continue;end
			intSurrY = vecSurrY(intSurrPix)+intY;
			if intSurrY < 1 || intSurrY > intSubY,continue;end
			
			%build pixel-wise distance matrix from correlations
			matMap = squeeze(matPixelCorrMat(intSurrY,intSurrX,:,:))-matMean;
			matCorr(intSurrPix,:) = circshift(matMap(:),[intSurrPix+floor(intSurrPixTot/2) 0]);
		end
		
		%transform correlation to distance
		matCorr(diag(diag(true(size(matCorr))))) = 0;
		matDist = 1-matCorr;
		matDist(isnan(matDist)) = 0;
		
		%cluster pixels in this neighborhood
		intMaxClustNum = 11;
		try
			[intNumberOfClusters,vecSilhouetteDistances] = doFastClustering(matDist,intMaxClustNum);
		catch
			intNumberOfClusters = 0;
		end
		if ~isempty(intNumberOfClusters) && intNumberOfClusters > 1 && intNumberOfClusters < intMaxClustNum
			matTree = linkage(matDist,'ward');
			vecClusterID = cluster(matTree,'maxclust',intNumberOfClusters);
			
			%build cluster picture
			%matGrid = getFillGrid(zeros(intSpatWin,intSpatWin),vecSurrY+intSubSW+1,vecSurrX+intSubSW+1,vecClusterID);
			%imagesc(matGrid);drawnow;
			
			%get this cluster and assign to cell array for all pixels that share the cluster
			indThisCluster = vecClusterID==vecClusterID(intThisPixIdx);
			vecClusterX = vecSurrX(indThisCluster)+intX;
			vecClusterY = vecSurrY(indThisCluster)+intY;
			vecOtherPixIdx = nan(1,numel(vecClusterX));
			for intOtherPix=1:numel(vecClusterX)
				vecOtherPixIdx(intOtherPix) = find(vecPixX==vecClusterX(intOtherPix) & vecPixY==vecClusterY(intOtherPix));
			end
			cellClust{intPixelIdx} = vecOtherPixIdx;
		end
	end
	%% build consensus model
	matBelongsToSameAsNeighbor = zeros(intSubY,intSubX);
	indClusterPixels = ~cellfun(@isempty,cellClust);
	vecClusterPixels = find(indClusterPixels);
	intClusterPixels = numel(vecClusterPixels);
	for intClustPixIdx = 1:intClusterPixels
		%get pixels for this cluster
		intClustPix = vecClusterPixels(intClustPixIdx);
		vecThisCluster = cellClust{intClustPix};
		
		%build map
		matGrid = getFillGrid(zeros(intSubY,intSubX,'logical'),vecPixY(vecThisCluster),vecPixX(vecThisCluster),true(numel(vecThisCluster),1));
		matPerimeter = bwperim(matGrid,8);
		matBelongsToSameAsNeighbor = matBelongsToSameAsNeighbor + matGrid-matPerimeter;
	end
	
	%%
	%build filter for cell body boundaries
	[matFiltX,matFiltY] = meshgrid(-3.5:3.5,-3.5:3.5);
	matGridXY = cat(3,matFiltX,matFiltY);
	
	%build center dip
	vecParams1 = [0 1 0 1 0 1 0];%[z-const, Amp, x0, wx, y0, wy, phi]
	matZ1 = get2DGauss(vecParams1,matGridXY);
	matZ1 = matZ1 ./ sum(matZ1(:));
	%build surround blob
	vecParams2 = [0 1 0 3 0 3 0];%[z-const, Amp, x0, wx, y0, wy, phi]
	matZ2 = get2DGauss(vecParams2,matGridXY);
	matZ2 = matZ2 ./ sum(matZ2(:));
	%combine
	matFilter = matZ2 - matZ1;
	
	%get initial locations for ring-detection algorithm
	matEdges = imfilter(matBelongsToSameAsNeighbor,matFilter, 'replicate');
	boolEdges = imbinarize(imnorm(matEdges),'adaptive','ForegroundPolarity','bright','Sensitivity',0.3);
	boolEdges = bwareaopen(boolEdges,6); %remove small blobs
	
	matFilled = imfill(matEdges,4);
	matCenters = matFilled - matEdges;
	
	
	%% find circles in proc im
	cellMethod = {'PhaseCode','TwoStage'};
	cellPolarity = {'bright','dark'};
	mapColor = parula(256);
	vecSmallRange = 2:9;
	matCentroidsAvgIm = [];
	for intRadIdx=1:numel(vecSmallRange)
		intSmallRad = vecSmallRange(intRadIdx);
		intLargeRad = round(intSmallRad*2);
		[centers,radii,metric] =  imfindcircles(imProc,[intSmallRad intLargeRad],...
			'ObjectPolarity',cellPolarity{1},'Method',cellMethod{2},'Sensitivity',0.8,'EdgeThreshold',0.1);
		vecColor = ceil(metric*size(mapColor,1));
		
		%assign to center list
		matCentroidsAvgIm = cat(1,matCentroidsAvgIm,centers);
	end
	
	
	% find circles in correlated activation
	imSegment = imnorm(matCenters);
	cellMethod = {'PhaseCode','TwoStage'};
	cellPolarity = {'bright','dark'};
	mapColor = parula(256);
	vecSmallRange = 2:9;
	matCentroidsCorrMat = [];
	for intRadIdx=1:numel(vecSmallRange)
		intSmallRad = vecSmallRange(intRadIdx);
		intLargeRad = round(intSmallRad*2);
		[centers,radii,metric] =  imfindcircles(imSegment,[intSmallRad intLargeRad],...
			'ObjectPolarity',cellPolarity{1},'Method',cellMethod{2},'Sensitivity',0.8,'EdgeThreshold',0.1);
		vecColor = ceil(metric*size(mapColor,1));
		
		%assign to center list
		matCentroidsCorrMat = cat(1,matCentroidsCorrMat,centers);
		
	end
	%% remove overlapping centroids
	matAllCentroids = cat(1,matCentroidsCorrMat,matCentroidsAvgIm);
	dblMergeThreshold = 3; %in subsampled pixels
	intCentroidNum = size(matAllCentroids,1);
	indKeepCentroids = true(intCentroidNum,1);
	intPointer = 1;
	while intPointer < intCentroidNum
		%get to-be-processed entries
		indSelectCentroids = indKeepCentroids;
		indSelectCentroids(1:intPointer) = false;
		
		%get this centroid's location
		vecLoc=matAllCentroids(intPointer,:);
		
		%get distance to other centroids
		vecDist = sqrt((vecLoc(1)-matAllCentroids(:,1)).^2 + (vecLoc(2)-matAllCentroids(:,2)).^2);
		
		%remove centroids that are too close
		indKeepCentroids(indSelectCentroids&(vecDist<dblMergeThreshold)) = false;
		intPointer = find(indKeepCentroids&indSelectCentroids,1);
	end
	
	%remove centroids
	matAllCentroids(~indKeepCentroids,:) = [];
	if boolPlot
	subplot(2,2,4)
	imshow(imProc);hold on;
	scatter(matAllCentroids(:,1),matAllCentroids(:,2),'bx');hold off;
	end
	%% loop through centroids and define boundaries
	%figure;
	intROIs = size(matAllCentroids,1);
	matMasksCorr = zeros(intSubY,intSubX,intROIs,'single');
	vecR2Corr = zeros(1,intROIs);
	matMasksIm = zeros(intSubY,intSubX,intROIs,'single');
	vecR2Im = zeros(1,intROIs);
	matFitParams = nan(7,intROIs);
	for intROI=1:intROIs
		%clf;
		%% get ROI x-y
		vecLoc = matAllCentroids(intROI,:);
		dblX = round(vecLoc(1));
		dblY = round(vecLoc(2));
		
		
		for intType = [1 2]
			try
				if boolPlot
					subplot(2,2,1);
					imshow(imProc);hold on;
					scatter(dblX,dblY,'rx');
				end
				
				%get location indices
				vecSelectX = dblX+[-intSubSW:intSubSW];
				vecSelectX(vecSelectX < 1 | vecSelectX > intSubX) = [];
				vecSelectY = dblY+[-intSubSW:intSubSW];
				vecSelectY(vecSelectY < 1 | vecSelectY > intSubY) = [];
				
				if intType == 1
					%corr-based
					%title('Corr')
					matMap = double(squeeze(matPixelCorrMat(dblY,dblX,:,:))-matMean);
					matMap(intSubSW+1,intSubSW+1) = matMap(intSubSW,intSubSW);
				else
					%im-based
					%title('Im')
					matMap = imProc(vecSelectY,vecSelectX);
				end
				if boolPlot
					subplot(2,2,2);
					imagesc(matMap);
				end
				%% calculate x-distro & y-distro as first guess for parameters
				%vecFitValsX = dblGainX*(normpdf(vecGrid,dblMuX,dblSigmaX)/dblConstantX); %proper function
				%vecFitVals = coeffvals(1)*exp(-((vecGrid-coeffvals(2))/coeffvals(3)).^2); %stupid function used by fit()
				%plot(fFit,vecGrid,vecBaseRemX');hold on;plot(vecGrid,vecFitVals,'bx');hold off;
				
				vecX = xmean(matMap,1);
				vecY = xmean(matMap,2);
				dblMin = max(min(vecX),min(vecY));
				vecBaseRemX = vecX-dblMin;
				vecBaseRemY = vecY-dblMin;
				
				%get x
				vecGridX = 1:size(matMap,2);
				fFitX = fit(vecGridX',vecBaseRemX','gauss1');
				vecCoeffsX = coeffvalues(fFitX);
				
				dblGainX = vecCoeffsX(1);
				dblMuX = ceil(numel(vecSelectX)/2);
				dblSigmaX = vecCoeffsX(3)/sqrt(2);
				dblConstantX = (1/sqrt(2*pi*dblSigmaX^2));
				
				% get y
				vecGridY = 1:size(matMap,1);
				fFitY = fit(vecGridY',vecBaseRemY,'gauss1');
				vecCoeffsY = coeffvalues(fFitY);
				
				dblGainY = vecCoeffsY(1);
				dblMuY = ceil(numel(vecSelectY)/2);
				dblSigmaY = vecCoeffsY(3)/sqrt(2);
				dblConstantY = (1/sqrt(2*pi*dblSigmaY^2));
				
				%check initial parameters
				if dblMuX < 1 || dblMuX > intSubSW*2,dblMuX = intSubSW;end
				if dblMuY < 1 || dblMuY > intSubSW*2,dblMuY = intSubSW;end
				
				%% fit 2D Gauss
				vecParam0 = [dblMin (dblGainX+dblGainY)/2 dblMuX dblSigmaX dblMuY dblSigmaY 0];%[z0, gain, x0, x-width, y0, y-width, angle(in rad)]
				[vecParams,matFit,resnorm,residual,exitflag] = getFitGauss2D(matMap,vecParam0);
				matFitParams(:,intROI) = vecParams;
				matFit = single(matFit-min(matFit(:)));
				if boolPlot
					subplot(2,2,3);
					imagesc(matFit);
				end
				%transform to original coordinates
				matThisMask = zeros(intSubY,intSubX,'single');
				matThisMask(vecSelectY,vecSelectX) = matFit;
				
				%calculate residuals
				matTot = matMap - mean(matMap(:));
				dblSS_Tot = sum(sum(matTot.^2));
				dblSS_Res = sum(residual.^2);
				dblR2 = 1 - (dblSS_Res / dblSS_Tot);
				if sqrt((vecParams(3)-dblMuX).^2 + (vecParams(5)-dblMuY).^2) > 5
					dblR2 = 0;
				end
				if boolPlot
					subplot(2,2,4);
					imagesc(matThisMask);
					title(sprintf('ROI %d/%d;R^2=%.3f',intROI,intROIs,dblR2));
					drawnow;
					pause(0.5);
				end
				
				
				%% assign to output
				if intType == 1
					%corr-based
					vecR2Corr(intROI) = dblR2;
					matMasksCorr(:,:,intROI) = matThisMask;
				else
					%im-based
					vecR2Im(intROI) = dblR2;
					matMasksIm(:,:,intROI) = matThisMask;
				end
			catch
			end
		end
	end
	vecR2Im(vecR2Im<0 | vecR2Im>1)=0;
	vecR2Corr(vecR2Corr<0 | vecR2Corr>1)=0;
	
	%% remove bad clusters, and remove overlaps
	matAutoCellDetectMasks = false(intSubY,intSubX,intROIs);
	
	%make mask to select outer rim
	matBaseMask = false(intSpatWin,intSpatWin);
	matBaseMask([1 2 end-1 end],:) = true;
	matBaseMask(:,[1 2 end-1 end]) = true;
	
	%sort clusters from high to low R2 for corr
	[dummy,vecReorder] = sort(vecR2Corr,'descend');
	for intSortedROI=1:numel(vecReorder)
		% get ROI idx and R2
		intOrigROI = vecReorder(intSortedROI);
		dblR2Corr = vecR2Corr(intOrigROI);
		dblR2Im = vecR2Im(intOrigROI);
		
		% get ROI x-y
		vecLoc = matAllCentroids(intOrigROI,:);
		dblX = round(vecLoc(1));
		dblY = round(vecLoc(2));
		
		%get location indices
		vecSelectX = dblX+[-intSubSW:intSubSW];
		vecSelectX(vecSelectX < 1 | vecSelectX > intSubX) = [];
		vecSelectY = dblY+[-intSubSW:intSubSW];
		vecSelectY(vecSelectY < 1 | vecSelectY > intSubY) = [];
		
		if dblR2Corr == 0 && dblR2Im == 0
			%% skip
			continue;
		elseif dblR2Corr >=  dblR2Im
			%% get data
			matMap = double(squeeze(matPixelCorrMat(dblY,dblX,:,:))-matMean);
			matMap(intSubSW+1,intSubSW+1) = matMap(intSubSW,intSubSW);
			matMask = matMasksCorr(vecSelectY,vecSelectX,intOrigROI);
			
			%% calculate optimal threshold with signal detection if correlation-based
			vecNoiseVals = matMap(matBaseMask);
			vecSignalVals = matMap(matMask>(mean(matMask(:))+std(matMask(:))));
			
			dblMu1 = mean(vecNoiseVals);
			dblSd1 = std(vecNoiseVals);
			
			dblMu2 = mean(vecSignalVals);
			dblSd2 = std(vecSignalVals);
			
			%build 1 - accuracy function
			fNegAcc = @(x) 1-((normcdf(x,dblMu1,dblSd1)+1-normcdf(x,dblMu2,dblSd2))/2);
			
			%get optimal criterion (highest accuracy)
			dblCriterion = fminsearch(fNegAcc,(dblMu1+dblMu2)/2);
			
			%get pixels
			matOptimMask = matMap > dblCriterion;
			
			%remove specks
			matOptimMask = bwareaopen(matOptimMask,32);
			
			%smooth
			matOptimMask = conv2(double(matOptimMask),[0 1 0; 1 1 1; 0 1 0],'same')>2;
			
			%fill
			matOptimMask = imfill(matOptimMask,'holes');
			
		elseif dblR2Corr <  dblR2Im
			%% get data
			matMap = imProc(vecSelectY,vecSelectX);
			matMask = matMasksIm(vecSelectY,vecSelectX,intOrigROI);
			
			%% use fit for mask
			matOptimMask = matMask;
			vecParams = matFitParams(:,intOrigROI);
			matOptimMask = (matOptimMask-min(matOptimMask(:)))/vecParams(2);
			matOptimMask = matOptimMask>(normpdf(1,0,1)/normpdf(0,0,1));
		end
		%% general processing
		%remove specks
		matOptimMask = bwareaopen(matOptimMask,32);
		
		%smooth
		matOptimMask = conv2(double(matOptimMask),[0 1 0; 1 1 1; 0 1 0],'same')>2;
		
		%fill
		matOptimMask = imfill(matOptimMask,'holes');
		
		%perform checks
		%check total size
		if sum(matOptimMask(:)) > (numel(matOptimMask)/2)
			continue;
		end
		%check number of objects
		[dummy,dummy2,intN] = bwboundaries(matOptimMask);
		if intN > 1
			continue;
		end
		
		%% check for overlap
		%transform to original coordinates
		matPaddedMask = false(intSubY,intSubX);
		matPaddedMask(vecSelectY,vecSelectX) = matOptimMask;
		
		matOverlap = bsxfun(@times,matAutoCellDetectMasks,matPaddedMask);
		if any(matOverlap(:))
			%build map
			matOverlapMap = any(matOverlap,3);
			vecOverlapPixels = squeeze(sum(sum(matOverlap,1),2));
			if max(vecOverlapPixels) > sum(matPaddedMask(:))/2
				%remove ROI if overlap is more than half
			else
				%remove overlapping pixels
				matPaddedMask = matPaddedMask - matOverlapMap;
				
				%assign
				matAutoCellDetectMasks(:,:,intOrigROI) = matPaddedMask;
			end
		else
			%assign
			matAutoCellDetectMasks(:,:,intOrigROI) = matPaddedMask;
		end
		
		%imshow(sum(matAutoCellDetectMasks,3));drawnow;
		%pause
	end
	
	%% remove bad ROIs
	indKeep = squeeze(sum(sum(matAutoCellDetectMasks,1),2))>0;
	matAutoCellDetectMasks(:,:,~indKeep) = [];
	
end
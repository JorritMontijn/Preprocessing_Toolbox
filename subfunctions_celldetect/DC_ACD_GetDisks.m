function [matMergedCentroids,vecMergedRadii] = DC_ACD_GetDisks(matImage,dblMergeThreshold,dblSensitivity,dblEdgeThreshold)
	
	%% input
	if ~exist('dblSensitivity','var') || isempty(dblSensitivity)
		dblSensitivity = 0.7;
	end
	if ~exist('dblEdgeThreshold','var') || isempty(dblEdgeThreshold)
		dblEdgeThreshold = 0.05;
	end
	
	%% detect
	warning('off','images:imfindcircles:warnForSmallRadius');
	cellMethod = {'PhaseCode','TwoStage'};
	cellPolarity = {'bright','dark'};
	vecSmallRange = 4:18;
	matCentroids = [];
	vecRadii = [];
	vecMetric = [];
	for intRadIdx=1:numel(vecSmallRange)
		intSmallRad = vecSmallRange(intRadIdx);
		intLargeRad = round(intSmallRad*2);
		[centers,radii,metric] =  imfindcircles(matImage,[intSmallRad intLargeRad],...
			'ObjectPolarity',cellPolarity{2},'Method',cellMethod{2},'Sensitivity',dblSensitivity,'EdgeThreshold',dblEdgeThreshold);
		
		%assign to list
		matCentroids = cat(1,matCentroids,centers);
		vecRadii = cat(1,vecRadii,radii);
		vecMetric = cat(1,vecMetric,metric);
	end
	warning('on','images:imfindcircles:warnForSmallRadius');
	
	%% decimate
	%pre-allocate merged list
	matMergedCentroids = nan(size(matCentroids));
	vecMergedRadii = nan(size(vecRadii));
	intMergeCounter = 0;
	
	%loop
	if ~exist('dblMergeThreshold','var') || isempty(dblMergeThreshold),dblMergeThreshold=10;end
	intCentroidNum = size(matCentroids,1);
	indKeepCentroids = true(intCentroidNum,1);
	intPointer = 1;
	while intPointer < intCentroidNum
		%get to-be-processed entries
		indSelectCentroids = indKeepCentroids;
		indSelectCentroids(1:(intPointer-1)) = false;
		
		%get this centroid's location
		vecLoc=matCentroids(intPointer,:);
		
		%get distance to other centroids
		vecDist = sqrt((vecLoc(1)-matCentroids(:,1)).^2 + (vecLoc(2)-matCentroids(:,2)).^2);
		
		%merge centroids that are too close
		indMergeCentroids = indSelectCentroids&(vecDist<dblMergeThreshold);
		indKeepCentroids(indMergeCentroids) = false;
		intPointer = find(indKeepCentroids&indSelectCentroids,1);
		%get constituent detections
		matTempCentroids =  matCentroids(indMergeCentroids,:);
		matTempRadii = vecRadii(indMergeCentroids,:);
		vecTempMetric = vecMetric(indMergeCentroids);
		%merge
		vecMergedCentroid = sum(bsxfun(@times,matTempCentroids,vecTempMetric),1)./sum(vecTempMetric);
		dblMergedRadius = sum(bsxfun(@times,matTempRadii,vecTempMetric),1)./sum(vecTempMetric);
		
		%assign
		intMergeCounter = intMergeCounter + 1;
		matMergedCentroids(intMergeCounter,:) = vecMergedCentroid;
		vecMergedRadii(intMergeCounter) = dblMergedRadius;
	end
	matMergedCentroids(intMergeCounter+1:end,:) = [];
	vecMergedRadii(intMergeCounter+1:end) = [];
	
end
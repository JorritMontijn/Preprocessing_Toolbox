function [vecDeleteObjects,matCannotBeFixed] = DC_ACD_RefineBorders(indForceAlgos)
	%refines cell boundaries
	
	%globals
	global sFig;
	global sRec;
	global sDC;
	
	
	%get algorithm
	strBoundaryType = sDC.metaData.cellBoundaryType{get(sFig.ptrListBoundaryType, 'Value')};
	if contains(strBoundaryType,'OGB','IgnoreCase',true)
		intDetectionType = 1;
	elseif contains(strBoundaryType,'GCaMP','IgnoreCase',true)
		intDetectionType = 2;
	elseif contains(strBoundaryType,'HoughFlood','IgnoreCase',true)
		intDetectionType = 3;
	end
	
	%get delete/keep switch
	strDeleteType = sDC.metaData.cellAllowDeletion{get(sFig.ptrListAllowDeletion, 'Value')};
	if contains(strDeleteType,'Delete','IgnoreCase',true)
		intDeleteType = 1;
	elseif contains(strDeleteType,'Keep','IgnoreCase',true)
		intDeleteType = 0;
	end
	
	%get refine type
	indUseRefineAlgos = false(1,2);
	strRefineType = sDC.metaData.cellRefineType{get(sFig.ptrListRefineType, 'Value')};
	if contains(strRefineType,'MultiTier+','IgnoreCase',true)
		indUseRefineAlgos(:) = true;
	elseif contains(strRefineType,'Corr','IgnoreCase',true)
		indUseRefineAlgos(1) = true;
	elseif contains(strRefineType,'Lum','IgnoreCase',true)
		indUseRefineAlgos(2) = true;
	end
	if exist('indForceAlgos','var')
		indUseRefineAlgos = indForceAlgos;
	end
	
	%msg
	clear cellText;
	cellText{1} = sprintf('Starting border refinement...');
	DC_updateTextInformation(cellText);
	
	%get all masks
	intObjects = numel(sDC.ROI);
	matMasks = cell2mat(reshape(arrayfun(@(x) logical(x.matMask),sDC.ROI,'UniformOutput',false),[1 1 intObjects]));
	matCannotBeFixed = [];
	
	%loop through ROIs and calculate overlap
	for intObject=1:intObjects
		%loop through ROIs and calculate overlap
		cellText{1} = sprintf('Refining border for object %d',intObject);
		DC_updateTextInformation(cellText);
	
		%get mask and compare with others
		matThisMask = sDC.ROI(intObject).matMask;
		
		%get overlaps
		matOverlaps = bsxfun(@and,matThisMask,matMasks);
		vecOverlaps = find(any(any(matOverlaps,1),2));
		%remove self
		vecOverlaps(vecOverlaps==intObject) = [];
		if isempty(vecOverlaps),continue;end
		
		for intOverlapIdx=1:numel(vecOverlaps)
			%keep track of fix
			boolFixed = false;
			intOtherObject = vecOverlaps(intOverlapIdx);
			
			%use pixel-wise correlation
			if indUseRefineAlgos(1) && isfield(sRec,'sPixResp') && isfield(sRec.sPixResp,'matPixelCorrMat')
				%get masks for overlap and both areas
				matThisMask = sDC.ROI(intObject).matMask;
				matOtherMask = sDC.ROI(intOtherObject).matMask;
				matThisOverlap = matThisMask & matOtherMask;
				matOnlyObject1 = matThisMask & ~matThisOverlap;
				matOnlyObject2 = matOtherMask & ~matThisOverlap;
				
				%check if one falls entirely in the other
				if sum(matOnlyObject1(:)) == 0 || sum(matOnlyObject2(:)) == 0
					%check if switch is set to keep
					if intDeleteType ~= 0
						%get internal consistency for object 1
						[dblMeanCorr1,dblSDCorr1] = DC_ACD_GetCorrPixForMask(matThisMask,sRec.sPixResp.matPixelCorrMat);
						
						%get internal consistency for object 2
						[dblMeanCorr2,dblSDCorr2] = DC_ACD_GetCorrPixForMask(matOtherMask,sRec.sPixResp.matPixelCorrMat);
						
						%get d'
						dblDprime = abs((dblMeanCorr1-dblMeanCorr2)/sqrt( ( (dblSDCorr1)^2 + (dblSDCorr2)^2)/2));
						
						%remove object with lowest consistency
						if dblDprime > 1 && dblMeanCorr1 > dblMeanCorr2 %delete 2
							DC_UpdateMask(intOtherObject,false(size(sDC.ROI(intOtherObject).matMask)));
							boolFixed = true;
						elseif dblDprime > 1 && dblMeanCorr1 < dblMeanCorr2 %delete 1
							DC_UpdateMask(intObject,false(size(sDC.ROI(intObject).matMask)));
							boolFixed = true;
						end
					end
				else
					%otherwise, check overlap and assign to 1 or 2
					[matMeanCorrOverlap1,matSDCorrOverlap1] = DC_ACD_GetCrossCorrPix(sRec.sPixResp.matPixelCorrMat,matThisOverlap,matThisMask);
					[matMeanCorrOverlap2,matSDCorrOverlap2] = DC_ACD_GetCrossCorrPix(sRec.sPixResp.matPixelCorrMat,matThisOverlap,matOtherMask);
					%get d'
					matDprime = abs((matMeanCorrOverlap1-matMeanCorrOverlap2)./sqrt( ( (matSDCorrOverlap1).^2 + (matSDCorrOverlap2).^2)/2));
					%calc regions	
					matRemoveFrom1 = (matMeanCorrOverlap2 > matMeanCorrOverlap1 & matDprime > 1) | matMeanCorrOverlap2 == matMeanCorrOverlap1;
					matRemoveFrom2 = (matMeanCorrOverlap1 > matMeanCorrOverlap2 & matDprime > 1) | matMeanCorrOverlap2 == matMeanCorrOverlap1;
					
					%remove from 1
					matMask = sDC.ROI(intObject).matMask;
					matMask(matRemoveFrom1) = false;
					DC_UpdateMask(intObject,matMask);
					
					%remove region 2
					matOtherMask = sDC.ROI(intOtherObject).matMask;
					matOtherMask(matRemoveFrom2) = false;
					DC_UpdateMask(intOtherObject,matOtherMask);
					
					if sum(matDprime(:) > 1) == sum(matThisOverlap(:))
						boolFixed = true;
					end
				end
			end
			if indUseRefineAlgos(2) %use luminance
				%get goodness of fit for both areas
				matThisMask = sDC.ROI(intObject).matMask;
				matOtherMask = sDC.ROI(intOtherObject).matMask;
				matThisOverlap = matThisMask & matOtherMask;
				matOnlyObject1 = matThisMask & ~matThisOverlap;
				matOnlyObject2 = matOtherMask & ~matThisOverlap;
				%get image
				matIm = sFig.cellIm{1}(:,:,2);
				%get luminance object 1
				vecPixLum1 = matIm(matOnlyObject1);
				vecPixLum2 = matIm(matOnlyObject2);
				dblLumMu1 = mean(vecPixLum1);
				dblLumMu2 = mean(vecPixLum2);
				dblLumSD1 = std(vecPixLum1);
				dblLumSD2 = std(vecPixLum2);
				%get d'
				dblDprime = abs((dblLumMu1-dblLumMu2)/sqrt( ( (dblLumSD1)^2 + (dblLumSD2)^2)/2));
				if abs(dblDprime) > 1
					if dblLumMu1 > dblLumMu2
						%remove from 1
						matMask = sDC.ROI(intOtherObject).matMask;
						matMask(matThisOverlap) = false;
						DC_UpdateMask(intOtherObject,matMask);
					else
						%remove from 1
						matMask = sDC.ROI(intObject).matMask;
						matMask(matThisOverlap) = false;
						DC_UpdateMask(intObject,matMask);
					end
					boolFixed = true;
				end
			end
			if ~boolFixed
				%add to unfixed list
				matCannotBeFixed(end+1,:) = [intObject intOtherObject]; %#ok<AGROW>
			end
		end
	end
	
	%recalculate borders if MultiTier+ is active
	vecDeleteObjects = [];
	if all(indUseRefineAlgos(:))
		matNewMasks = cell2mat(reshape(arrayfun(@(x) logical(x.matMask),sDC.ROI,'UniformOutput',false),[1 1 numel(sDC.ROI)]));
		vecObjectSizes = squeeze(sum(sum(matNewMasks,2),1));
		dblPixNumThresh = 30;
		vecCheckObjects = find(vecObjectSizes<dblPixNumThresh);
		if intDetectionType == 3
			for intObjectIdx=1:numel(vecCheckObjects)
				intObject = vecCheckObjects(intObjectIdx);
				%redetect borders
				DC_UpdateMask(intObject,DC_ACD_HoughFloodPlus(sDC.ROI(intObject).intCenterX, sDC.ROI(intObject).intCenterY));
			end
		else
			[intDetectObjects,indDelete] = DC_detectObjects(vecCheckObjects,true);
		end
		
		%rerun same function, but with no redetection
		indForceAlgos = false(size(indUseRefineAlgos));
		indForceAlgos(1) = true;
		[vecDeleteObjects2,matCannotBeFixed] = DC_ACD_RefineBorders(indForceAlgos);
		vecDeleteObjects = cat(2,vecDeleteObjects,vecDeleteObjects2);
	end
	
	%final check
	matNewMasks = cell2mat(reshape(arrayfun(@(x) logical(x.matMask),sDC.ROI,'UniformOutput',false),[1 1 numel(sDC.ROI)]));
	vecObjectSizes = squeeze(sum(sum(matNewMasks,2),1));
	dblPixNumThresh = 10;
	vecCheckObjects = find(vecObjectSizes<dblPixNumThresh);
	if intDeleteType == 0
		for intObjectIdx=1:numel(vecCheckObjects)
			intObject = vecCheckObjects(intObjectIdx);
			matMask = false(size(sDC.ROI(intObject).matMask));
			matMask([-1:1]+round(sDC.ROI(intObject).intCenterY),[-1:1]+round(sDC.ROI(intObject).intCenterX)) = true;
			DC_UpdateMask(intObject,matMask);
		end
	else
		vecDeleteObjects = cat(2,vecDeleteObjects,vecCheckObjects(:)');
	end
	
	%delete objects
	vecDeleteObjects = unique(vecDeleteObjects);
	if numel(vecDeleteObjects) > 0
		DC_removeObject(vecDeleteObjects)
	end
end
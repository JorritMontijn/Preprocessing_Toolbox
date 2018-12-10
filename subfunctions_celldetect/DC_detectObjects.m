function [intDetectObjects,indDelete] = DC_detectObjects(vecObjects,boolForce)
	%detects cell boundaries
	
	%globals
	global sFig;
	global sRec;
	global sDC;
	
	%check inputs
	if exist('vecObjects','var'),vecObjects = vecObjects(:)';end
	if ~exist('vecObjects','var') || isempty(vecObjects),vecObjects = 1:numel(sDC.ROI);end
	if ~exist('boolForce','var'),boolForce=false;end
	
	%get general size
	dblPixelSize = sRec.xml.sData.dblActualVoxelSizeX / 1000; % micrometers
	
	%calc number of objects
	intObjects = numel(sDC.ROI);
	strBoundaryType = sDC.metaData.cellBoundaryType{get(sFig.ptrListBoundaryType, 'Value')};
	if contains(strBoundaryType,'OGB','IgnoreCase',true)
		intDetectionType = 1;
		boolInvertNeurons = false;
	elseif contains(strBoundaryType,'GCaMP','IgnoreCase',true)
		intDetectionType = 2;
		boolInvertNeurons = true;
	elseif contains(strBoundaryType,'HoughFlood','IgnoreCase',true)
		intDetectionType = 3;
		boolInvertNeurons = true;
	end
	
	%get refine type
	strRefineType = sDC.metaData.cellRefineType{get(sFig.ptrListRefineType, 'Value')};
	if contains(strRefineType,'MultiTier+','IgnoreCase',true)
		intUseAlgo = 1;
	elseif contains(strRefineType,'Corr','IgnoreCase',true)
		intUseAlgo = 1;
	elseif contains(strRefineType,'Lum','IgnoreCase',true)
		intUseAlgo = 2;
	end
	
	%get original image
	intSelectedImage = get(sFig.ptrListSelectImage,'Value');
	if intSelectedImage > 2
		imDetect = sFig.cellIm{intSelectedImage};
		imHough = mean(imDetect,3);
	else
		imHough = [];
		imDetect = sFig.cellIm{2};
	end
	imGreen = imDetect(:,:,2);
	imRed = imDetect(:,:,1);
	
	%change im1D if necessary
	if intUseAlgo == 2
		imUse = sFig.cellIm{1}(:,:,2);
		imUse = (imUse + circshift(imUse,[1 0]))/2;
		imHough = 1-imnorm(imUse);
	end
	
	%pre-allocate
	intDetectObjects = 0;
	vecOriginalIndex = nan(1,intObjects);
	for intObject = vecObjects
		if isempty(sDC.ROI(intObject).matMask) || boolForce
			%increment counter
			intDetectObjects = intDetectObjects + 1;
			vecOriginalIndex(intDetectObjects) = intObject;
		end
	end
	cellImages = cell(1,intDetectObjects);
	vecOriginalIndex = vecOriginalIndex(1:intDetectObjects);
	
	%set counters
	intDetectObject = 0;
	indDelete = []; %for backward compatibility
	
	% get cell mask
	for intObject = vecOriginalIndex
		%increment counter
		intDetectObject = intDetectObject + 1;
		
		%get object data
		intX = sDC.ROI(intObject).intCenterX;
		intY = sDC.ROI(intObject).intCenterY;
		
		%detect object
		if intDetectionType == 1 || intDetectionType == 2
			%% old detection
			%get type
			intType = sDC.ROI(intObject).intType;
			strType = sDC.metaData.cellType{intType};
			
			%assign color
			boolNeuron = false;
			dblROISize = sDC.metaData.dblExpectedCellSize;
			if strcmp(strType,'neuron')
				boolNeuron = true;
				boolInvert = boolInvertNeurons;
				dblThresholdPercentage = 0.6;
				imGaussDetect = imGreen;
			elseif strcmp(strType,'astrocyte')
				boolNeuron = true;
				boolInvert = false;
				dblThresholdPercentage = 0.65;
				imGaussDetect = imRed;
			elseif strcmp(strType,'bloodvessel')
				boolInvert = true;
				dblThresholdPercentage = 0.7;
				imGaussDetect = imGreen;
			elseif strcmp(strType,'neuropil')
				boolInvert = false;
				dblThresholdPercentage = 0.3;
				imGaussDetect = imGreen;
			elseif strcmp(strType,'PV')
				boolNeuron = true;
				boolInvert = boolInvertNeurons;
				dblThresholdPercentage = 0.55;
				imGaussDetect = imGreen;
			elseif strcmp(strType,'SOM')
				boolNeuron = true;
				boolInvert = boolInvertNeurons;
				dblThresholdPercentage = 0.6;
				imGaussDetect = imGreen;
			elseif strcmp(strType,'VIP')
				boolNeuron = true;
				boolInvert = boolInvertNeurons;
				dblThresholdPercentage = 0.6;
				imGaussDetect = imGreen;
			else
				boolInvert = false;
				dblThresholdPercentage = 0.7;
			end
			%vecColor = sDC.metaData.cellColor{intType};
			
			%get mask
			imDetect = DC_ACD_maskCellSurround(imGaussDetect, intX, intY, dblPixelSize, dblROISize, boolInvert);
			% get thresholded contour
			if intDetectionType == 2 && boolNeuron
				%flatten core
				imCenter = DC_ACD_thresholdContour(imDetect, dblThresholdPercentage);
				imGaussDetect(imCenter>0) = 1;
				
				%detect surrounding border
				imDetectBorder = DC_ACD_maskCellSurround(imGaussDetect, intX, intY, dblPixelSize, dblROISize*1.3, false);
				imBorder = DC_ACD_thresholdContour(imDetectBorder, 0.70);
				
				%dilate somewhat to include whole object
				cellImages{intDetectObject} = imdilate(imBorder,strel('disk', 3, 0),'same');
			else
				cellImages{intDetectObject} = DC_ACD_thresholdContour(imDetect, dblThresholdPercentage);
			end
			
		elseif intDetectionType == 3
			%HoughFlood+
			DC_updateTextInformation({sprintf('Processing object %d/%d...',intDetectObject,numel(vecOriginalIndex))});
			matMask = DC_ACD_HoughFloodPlus(intX, intY, imHough);
			
			%get centroid
			[vecYX] = calcCenterOfMass(matMask);
			
			%assign to structure
			sDC.ROI(intObject).intCenterX = vecYX(2);
			sDC.ROI(intObject).intCenterY = vecYX(1);
			sDC.ROI(intObject).matMask = matMask;
			sDC.ROI(intObject).dblRadius = [];
			
			%delete marker
			if sFig.sObject(intObject).drawn == 1 && isfield(sFig.sObject(intObject),'handles') && isfield(sFig.sObject(intObject).handles,'marker') && ~isempty(sFig.sObject(intObject).handles.marker)
				delete(sFig.sObject(intObject).handles.marker);
			end
			
			%set redraw flags
			sFig.sObject(intObject).drawn = 0;
			sFig.sObject(intObject).handles.marker = [];
			sFig.sObject(intObject).handles.lines = [];
		end
	end
	
	if intDetectionType < 3
		if intDetectObject == 0
			indDelete = [];
			return;
		end
		
		% remove overlap
		[ImRO, oList, indDelete] = DC_ACD_removeOverlap(cellImages);
		for intDetectObject = 1:intDetectObjects
			%get mask data
			[matMask,nPolygons] = bwlabel(ImRO{intDetectObject},4) ;  % also 8 possible
			
			%get properties of identified objects
			cellBasic  = regionprops(matMask, 'Basic') ;
			%cellPerimeter = regionprops(matMask, 'ConvexHull') ;
			intObject = vecOriginalIndex(intDetectObject);
			
			%assign to structure
			sDC.ROI(intObject).intCenterX = cellBasic(1).Centroid(1);
			sDC.ROI(intObject).intCenterY = cellBasic(1).Centroid(2);
			
			%sDC.ROI(intObject).matPerimeter = cellPerimeter(1).ConvexHull ;
			sDC.ROI(intObject).matMask = matMask;
			
			%delete marker
			if sFig.sObject(intObject).drawn == 1 && isfield(sFig.sObject(intObject),'handles') && isfield(sFig.sObject(intObject).handles,'marker') && ~isempty(sFig.sObject(intObject).handles.marker)
				delete(sFig.sObject(intObject).handles.marker);
			end
			
			%set redraw flags
			sFig.sObject(intObject).drawn = 0;
			sFig.sObject(intObject).handles.marker = [];
			sFig.sObject(intObject).handles.lines = [];
		end
	end
end
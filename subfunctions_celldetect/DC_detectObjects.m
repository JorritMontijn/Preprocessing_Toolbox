function [intDetectObjects,indDelete] = DC_detectObjects(vecObjects)
	%detects cell boundaries
	
	%globals
	global sFig;
	global sRec;
	global sDC;
	
	%get general size
	dblPixelSize = sRec.xml.sData.dblActualVoxelSizeX / 1000; % micrometers
	
	%calc number of objects
	intObjects = numel(sDC.ROI);
	strBoundaryType = sDC.metaData.cellBoundaryType{get(sFig.ptrListBoundaryType, 'Value')};
	if strcmpi(strBoundaryType,'OGB')
		intDetectType = 1;
		boolInvertNeurons = false;
	elseif strcmpi(strBoundaryType,'GCaMP')
		intDetectType = 2;
		boolInvertNeurons = true;
	end
	if ~exist('vecObjects','var') || isempty(vecObjects),vecObjects = 1:intObjects;end
	
	%get original image
	imDetect = sFig.cellIm{2};
	imGreen = imDetect(:,:,2);
	imRed = imDetect(:,:,1);
	
	%pre-allocate
	intDetectObjects = 0;
	vecOriginalIndex = nan(1,intObjects);
	for intObject = vecObjects
		if isempty(sDC.ROI(intObject).matMask)
			
			%increment counter
			intDetectObjects = intDetectObjects + 1;
			vecOriginalIndex(intDetectObjects) = intObject;
		end
	end
	cellImages = cell(1,intDetectObjects);
	vecOriginalIndex = vecOriginalIndex(1:intDetectObjects);
	
	%set counters
	intDetectObject = 0;
	
	% get cell mask
	for intObject = vecOriginalIndex
		%increment counter
		intDetectObject = intDetectObject + 1;
		
		%get object data
		intX = sDC.ROI(intObject).intCenterX;
		intY = sDC.ROI(intObject).intCenterY;
		
		%get type
		intType = sDC.ROI(intObject).intType;
		strType = sDC.metaData.cellType{intType};
		
		%assign color
		boolNeuron = false;
		dblROISize = sDC.metaData.dblExpectedCellSize;
		if strcmp(strType,'neuron')
			boolNeuron = true;
			boolInvert = boolInvertNeurons;
			dblThresholdPercentage = 0.7;
			im1D = imGreen;
		elseif strcmp(strType,'astrocyte')
			boolInvert = false;
			dblThresholdPercentage = 0.78;
			im1D = imRed;
		elseif strcmp(strType,'bloodvessel')
			boolInvert = true;
			dblThresholdPercentage = 0.8;
			im1D = imGreen;
		elseif strcmp(strType,'neuropil')
			boolInvert = false;
			dblThresholdPercentage = 0.3;
			im1D = imGreen;
		elseif strcmp(strType,'PV')
			boolNeuron = true;
			boolInvert = boolInvertNeurons;
			dblThresholdPercentage = 0.7;
			im1D = imGreen;
		elseif strcmp(strType,'SOM')
			boolNeuron = true;
			boolInvert = boolInvertNeurons;
			dblThresholdPercentage = 0.7;
			im1D = imGreen;
		elseif strcmp(strType,'VIP')
			boolNeuron = true;
			boolInvert = boolInvertNeurons;
			dblThresholdPercentage = 0.7;
			im1D = imGreen;
		else
			boolInvert = false;
			dblThresholdPercentage = 0.7;
		end
		%vecColor = sDC.metaData.cellColor{intType};
		
		%detect object
		imDetect = DC_ACD_maskCellSurround(im1D, intX, intY, dblPixelSize, dblROISize, boolInvert);
		
		% get thresholded contour
		if intDetectType == 2 && boolNeuron
			%flatten core
			imCenter = DC_ACD_thresholdContour(imDetect, dblThresholdPercentage);
			im1D(imCenter>0) = 1;
			
			%detect surrounding border
			imDetectBorder = DC_ACD_maskCellSurround(im1D, intX, intY, dblPixelSize, dblROISize*1.3, false);
			imBorder = DC_ACD_thresholdContour(imDetectBorder, 0.70);
			
			%dilate somewhat to include whole object
			cellImages{intDetectObject} = imdilate(imBorder,strel('disk', 3, 0),'same');
		else
			cellImages{intDetectObject} = DC_ACD_thresholdContour(imDetect, dblThresholdPercentage);
		end
	end
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
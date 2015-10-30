function sDC = doTransformCDGtosDC(CDG)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	vecSize = size(CDG.I{1});
	matMaskBase = false(vecSize(1),vecSize(2));
	
	sDC.strRecFile = '';
	sDC.strRecPath = '';
	
	sDC.metaData = struct;
	sDC.metaData.cellType = CDG.cellTypes;
	sDC.metaData.cellColor = {[0.5,0.1,0.5],[0.5,0.5,0],[1,0,0],[0.5,0,0.5],[0,0,1],[0,0.5,1],[0,1,1],[1,1,1],[1,1,1];};
	sDC.metaData.vecNeurons=[1 5 6 7];
	sDC.metaData.dblExpectedCellSize = 10;
	
	sDC.ROI = struct;
	for intROI=1:numel(CDG.cells)
		matMask = matMaskBase;
		for intPoint=1:size(CDG.cells(intROI).Body,1)
			int1=CDG.cells(intROI).Body(intPoint,1);
			int2=CDG.cells(intROI).Body(intPoint,2);
			matMask(int1,int2) = true;
		end
		sDC.ROI(intROI).intType = CDG.cells(intROI).type;
		sDC.ROI(intROI).intCenterX = mean([min(CDG.cells(intROI).Body(:,2)) max(CDG.cells(intROI).Body(:,2))]);
		sDC.ROI(intROI).intCenterY = mean([min(CDG.cells(intROI).Body(:,1)) max(CDG.cells(intROI).Body(:,1))]);
		sDC.ROI(intROI).matPerimeter = CDG.cells(intROI).Perimeter;
		sDC.ROI(intROI).matMask = matMask;
	end
end
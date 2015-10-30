function sDC = DC_populateStructure()
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%path locations
	sDC.metaData.strProcessedPath = 'D:\Data\Processed\imagingdata';
	sDC.metaData.strRawPath = 'D:\Data\Raw\imagingdata';
	
	%cell types:
	sDC.metaData.cellPixRespType{1} = 'None';
	sDC.metaData.cellPixRespType{2} = 'Selectivity';
	sDC.metaData.cellPixRespType{3} = 'Activation';
	
	sDC.metaData.cellType{1} = 'neuron';
	sDC.metaData.cellType{2} = 'astrocyte';
	sDC.metaData.cellType{3} = 'bloodvessel';
	sDC.metaData.cellType{4} = 'neuropil';
	sDC.metaData.cellType{5} = 'PV';
	sDC.metaData.cellType{6} = 'SOM';
	sDC.metaData.cellType{7} = 'VIP';
	sDC.metaData.cellType{8} = 'empty';
	sDC.metaData.cellType{9} = 'ToBeDeleted';
	
	sDC.metaData.cellPresence{1} = 'present';
	sDC.metaData.cellPresence{2} = 'absent';
	sDC.metaData.cellPresence{3} = 'filled';

	sDC.metaData.cellRespType{1} = 'tun+resp';
	sDC.metaData.cellRespType{2} = 'tuned';
	sDC.metaData.cellRespType{3} = 'responsive';
	sDC.metaData.cellRespType{4} = 'silent';
	
	sDC.metaData.cellDrawType{1} = 'Border';
	sDC.metaData.cellDrawType{2} = 'Centroid';
	
	sDC.metaData.cellBoundaryType{1} = 'OGB';
	sDC.metaData.cellBoundaryType{2} = 'GCaMP';
		
	sDC.metaData.cellColor{1} = [1.0 0.0 0.0]; %neuron
	sDC.metaData.cellColor{2} = [0.5 0.5 0.0]; %astrocyte
	sDC.metaData.cellColor{3} = [0.0 0.0 1.0]; %bloodvessel
	sDC.metaData.cellColor{4} = [0.0 1.0 0.0]; %neuropil
	sDC.metaData.cellColor{5} = [1.0 0.5 0.0]; %PV
	sDC.metaData.cellColor{6} = [1.0 0.5 0.0]; %SOM
	sDC.metaData.cellColor{7} = [1.0 0.5 0.0]; %VIP
	sDC.metaData.cellColor{8} = [0.0 0.0 0.0]; %other
	sDC.metaData.cellColor{9} = [0.0 0.0 0.0]; %other
	sDC.metaData.vecNeurons = [1 5 6 7]; %which ROIs are neurons
	sDC.metaData.dblExpectedCellSize = 10; %microns
	sDC.metaData.intROIDisplacementX = 0;
	sDC.metaData.intROIDisplacementY = 0;
	sDC.metaData.intROIAssignedDispX = 0;
	sDC.metaData.intROIAssignedDispY = 0;
	
	sDC.ROI = [];
	%sDC.ROI(intObject).matMask
	%sDC.ROI(intObject).intType
	%sDC.ROI(intObject).matPerimeter
	%sDC.ROI(intObject).intCenterX
	%sDC.ROI(intObject).intCenterY
end
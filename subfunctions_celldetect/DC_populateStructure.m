function sDC = DC_populateStructure(sDC)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%check metadata
	if ~isfield(sDC,'metaData')
		sDC.metaData = struct;
	end
	
	%for backward compatibility
	
	
	%path locations
	sDC.metaData.strProcessedPath = 'D:\Data\Processed\imagingdata';
	sDC.metaData.strRawPath = 'D:\Data\Raw\imagingdata';
	
	%pixel response maps
	if ~isfield(sDC.metaData,'cellPixRespType')
	sDC.metaData.cellPixRespType{1} = 'None';
	sDC.metaData.cellPixRespType{2} = 'Selectivity';
	sDC.metaData.cellPixRespType{3} = 'Activation';
	end
	
	%cell presence
	if ~isfield(sDC.metaData,'cellPresence')
		sDC.metaData.cellPresence{1} = 'present';
		sDC.metaData.cellPresence{2} = 'absent';
		sDC.metaData.cellPresence{3} = 'filled';
	end
	
	%responsiveness list
	if ~isfield(sDC.metaData,'cellRespType')
		sDC.metaData.cellRespType{1} = 'tun+resp';
		sDC.metaData.cellRespType{2} = 'tuned';
		sDC.metaData.cellRespType{3} = 'responsive';
		sDC.metaData.cellRespType{4} = 'silent';
	end
	
	%draw types
	if ~isfield(sDC.metaData,'cellDrawType')
		sDC.metaData.cellDrawType{1} = 'Border';
		sDC.metaData.cellDrawType{2} = 'Centroid';
	end
	
	%algorithm
	if ~isfield(sDC.metaData,'cellBoundaryType')
		sDC.metaData.cellBoundaryType{1} = 'OGB (legacy)';
		sDC.metaData.cellBoundaryType{2} = 'GCaMP (legacy)';
		sDC.metaData.cellBoundaryType{3} = 'HoughFlood+';
	end
	
	%overlap metric
	if ~isfield(sDC.metaData,'cellRefineType')
		sDC.metaData.cellRefineType{1} = 'MultiTier+';
		sDC.metaData.cellRefineType{2} = 'CorrPix';
		sDC.metaData.cellRefineType{3} = 'Luminance';
	end
	
	%allow deletions
	if ~isfield(sDC.metaData,'cellAllowDeletion')
		sDC.metaData.cellAllowDeletion{1} = 'Delete Obsolete';
		sDC.metaData.cellAllowDeletion{2} = 'Keep Obsolete';
	end
	
	%neuron types and draw colors
	if ~isfield(sDC.metaData,'cellType') ||  ~isfield(sDC.metaData,'cellColor') || ~isfield(sDC.metaData,'vecNeurons')
		sDC.metaData.cellType{1} = 'neuron';
		sDC.metaData.cellType{2} = 'astrocyte';
		sDC.metaData.cellType{3} = 'bloodvessel';
		sDC.metaData.cellType{4} = 'neuropil';
		sDC.metaData.cellType{5} = 'PV';
		sDC.metaData.cellType{6} = 'SOM';
		sDC.metaData.cellType{7} = 'VIP';
		sDC.metaData.cellType{8} = 'empty';
		sDC.metaData.cellType{9} = 'ToBeDeleted';
		
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
	end
	
	%displacement
	if ~isfield(sDC.metaData,'intROIDisplacementX')
		sDC.metaData.intROIDisplacementX = 0;
		sDC.metaData.intROIDisplacementY = 0;
		sDC.metaData.intROIAssignedDispX = 0;
		sDC.metaData.intROIAssignedDispY = 0;
	end
	
	%ROI structure
	if ~isfield(sDC,'ROI')
		sDC.ROI = [];
	end
	
end
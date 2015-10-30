function [MaskCell, MaskCellNeuropil] = getROImasksPrePro(sRec,sDC)


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% prepare script variables

	% number of rois
	numROIs = numel(sDC.ROI);
	
	% structure element 1 micron, 2 micron and 5 micron
	SEin = strel( 'disk', round(2 / (sRec.xml.sData.dblActualVoxelSizeX/1000)) ); %inner ring (ignored)
	SEout = strel( 'disk', round(5 / (sRec.xml.sData.dblActualVoxelSizeX/1000)) ); %outer ring (used for fluorescence correction)

	% mask variables
	MaskCell = cell(1,numROIs);
	MaskCellNeuropil = cell(1,numROIs);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% create masks for every roi
	
	for c = 1:numROIs

		%rename
		MaskCell{c} = logical(sDC.ROI(c).matMask);
	end
	
	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% create cell-removed mask

	% initialize mask variable
	MaskCellRemoved = true(sRec.sProcLib.y, sRec.sProcLib.x);

	% remove cells
	for c = 1:numROIs

		% get outer region of roi
		DilatedCellMask = imdilate( MaskCell{c}, SEin );

		% remove from masked area
		MaskCellRemoved( DilatedCellMask==1 ) = 0;
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% get neuropil area per roi

	% loop per roi
	for c = 1:numROIs

		% get outer region of roi
		InnerMask = imdilate( MaskCell{c}, SEin );

		% get outer region of surrounding neuropil
		MaskCellNeuropil{c} = imdilate( MaskCell{c}, SEout );

		% get area between inner and outer mask
		MaskCellNeuropil{c}( InnerMask==1 ) = 0;

		% remove area that is covered by other extended ROIs (ROIs+2micron)
		MaskCellNeuropil{c}( MaskCellRemoved==0 ) = 0;
	end
end


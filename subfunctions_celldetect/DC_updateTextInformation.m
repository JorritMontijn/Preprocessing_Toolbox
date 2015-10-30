function DC_updateTextInformation(varargin)
	%update cell information window
	global sFig;
	global sDC;
	
	%check if data has been loaded
	if isempty(sDC) || isempty(sFig)
		return;
	else
		try
			cellOldText = get(sFig.ptrTextInformation, 'string');
		catch
			return;
		end
	end
	
	%check if msg is supplied, otherwise display cell data
	intObject = 0;
	intSubTypeNr = 0;
	if nargin > 0
		cellText = varargin{1};
	else
		% if one cell is selected, display info
		if length(sFig.vecSelectedObjects) == 1
			intObject = sFig.vecSelectedObjects;
			intSubTypeNr = getSubTypeNr(sDC,intObject);
			if isfield(sDC.ROI(intObject),'intPresence') && isscalar(sDC.ROI(intObject).intPresence),intPresence=sDC.ROI(intObject).intPresence;else intPresence=1;end
			if isfield(sDC.ROI(intObject),'intRespType') && isscalar(sDC.ROI(intObject).intRespType),intRespType=sDC.ROI(intObject).intRespType;else intRespType=1;end
			set(sFig.ptrEditSelect,'string',num2str(intObject));
			set(sFig.ptrEditSelect,'string',num2str(intSubTypeNr));
			cellText{1} = ['object: ' num2str(intObject) ' ; subtype nr: ' num2str(intSubTypeNr)];
			cellText{2} = ['type: ' sDC.metaData.cellType{sDC.ROI(intObject).intType}];
			cellText{3} = ['presence: ' sDC.metaData.cellPresence{intPresence}];
			cellText{4} = ['responsiveness: ' sDC.metaData.cellRespType{intRespType}];
			cellText{5} = ['size: ' num2str(sum(sDC.ROI(intObject).matMask(:))) 'pixels'];
			cellText{6} = ['x-center: ' num2str(sDC.ROI(intObject).intCenterX)];
			cellText{7} = ['y-center: ' num2str(sDC.ROI(intObject).intCenterY)];
		elseif sum(sFig.vecSelectedObjects) > 1
			cellText{1} = ['number of selected objects: ' num2str(length(sFig.vecSelectedObjects))] ;
			cellText{2} = [''] ;
			cellText{3} = ['more than 1 cell selected, no info'] ;
		elseif sum(sFig.vecSelectedObjects) < 1
			cellText{1} = ['no cells selected'] ;
			cellText{2} = [''] ;
			cellText{3} = [''] ;
		end
		%set selection boxes
		set(sFig.ptrEditSelect,'string',num2str(intObject));
		set(sFig.ptrEditSubSelect,'string',num2str(intSubTypeNr));
	end
	
	
	
	%set figure text
	set(sFig.ptrTextInformation, 'string', cellText );
	drawnow;
end
function intSubTypeNr = getSubTypeNr(sDC,intObject)
	intType = sDC.ROI(intObject).intType;
	intSubTypeNr = 1;
	for intObjectCounter=1:(intObject-1)
		if intType == sDC.ROI(intObjectCounter).intType
			intSubTypeNr = intSubTypeNr + 1;
		end
	end
end
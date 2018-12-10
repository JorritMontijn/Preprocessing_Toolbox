function DC_updateTextInformation(varargin)
	%update cell information window
	global sFig;
	global sDC;
	
	%check if data has been loaded
	if isempty(sFig) || (isempty(sDC) && nargin == 0)
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
			cellText{1} = sprintf('Object: %d; subtype: %d',intObject,intSubTypeNr);
			cellText{2} = sprintf('Type: %s',sDC.metaData.cellType{sDC.ROI(intObject).intType});
			cellText{3} = sprintf('Presence: %s',sDC.metaData.cellPresence{intPresence});
			cellText{4} = sprintf('Responsiveness: %s',sDC.metaData.cellRespType{intRespType});
			cellText{5} = sprintf('Size: %d pixels',sum(sDC.ROI(intObject).matMask(:)));
			cellText{6} = sprintf('Center: x=%.2f, y=%.2f',sDC.ROI(intObject).intCenterX,sDC.ROI(intObject).intCenterY);
		elseif sum(sFig.vecSelectedObjects) > 1
			cellText{1} = sprintf('Number of selected objects: %d',length(sFig.vecSelectedObjects));
			cellText{2} = '';
			cellText{3} = 'More than 1 cell selected, no info';
		elseif sum(sFig.vecSelectedObjects) < 1
			cellText{1} = 'No cells selected';
		end
		%set selection boxes
		set(sFig.ptrEditSelect,'string',num2str(intObject));
		set(sFig.ptrEditSubSelect,'string',num2str(intSubTypeNr));
	end
	%{
	%add time stamp
	for intLine=1:numel(cellText)
		if ~isempty(cellText{intLine}) && ~strcmpi(cellText{intLine}(1),'[') 
			cellText{intLine} = [sprintf('[%s] ',getTime),cellText{intLine}];
		end
	end
	
	%set figure text
	cellTextOld = get(sFig.ptrTextInformation, 'string');
	cellText = cat(1,cellText,{''},cellTextOld);
	%}
	if numel(cellText) > 6,cellText(7:end) = [];end
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
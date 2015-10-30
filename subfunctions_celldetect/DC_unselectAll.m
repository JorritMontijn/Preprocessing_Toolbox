function boolChange = DC_unselectAll()
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
	%get globals
	global sDC;
	global sFig;
	global sRec;
	
	%check if data has been loaded
	if isempty(sDC) || isempty(sFig) || isempty(sRec)
		return;
	else
		try
			intType = get(sFig.ptrListCellType, 'Value');
			axesSize = get(sFig.ptrAxesHandle, 'Position');
			intPresence = get(sFig.ptrListPresence, 'Value');
		catch
			return;
		end
	end
	
	%loop through all selected objects to deselect them
	boolChange = false;
	for intObject = sFig.vecSelectedObjects
		if isfield(sDC.ROI,'matMask') && ~isempty(sDC.ROI(intObject).matMask) && sFig.sObject(intObject).drawn
			% remove old drawing and set drawn flag
			if isfield(sFig.sObject(intObject).handles,'lines')
				for p = 1:length(sFig.sObject(intObject).handles.lines)
					try delete(sFig.sObject(intObject).handles.lines(p));catch, end
				end
				sFig.sObject(intObject).handles.lines = [];
			else
				p = [];
			end
			if isempty(p) && isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX) && sFig.sObject(intObject).drawn == 1
				try delete(sFig.sObject(intObject).handles.marker);catch, end
			end
			sFig.sObject(intObject).drawn = 0;
			
			% end loop
			boolChange = true;
		end
	end
	
	%clear selection vector
	sFig.vecSelectedObjects = [];
end


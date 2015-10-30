function DC_removeObject(vecRemObjects)
	%removes objects from sDC structure
	
	%get globals
	global sDC;
	global sFig;
	
	%check if data has been loaded
	if isempty(sDC) || isempty(sFig)
		return;
	else
		try
			intImSelected = get(sFig.ptrListSelectImage,'Value');
		catch
			return;
		end
	end
	
	%get number of objects
	intObjects = numel(sDC.ROI);
	if intObjects == 0
		vecRemObjects = [];
	else
		%make vector of all objects
		vecAllObjects = 1:intObjects;
		
		%make index vector of which ones to remove
		vecRemIndex = ismember(vecAllObjects,vecRemObjects);
		
		for intObject=find(vecRemIndex)
			% remove old drawing and set drawn flag
			if isfield(sFig.sObject(intObject).handles,'lines')
				for p = 1:length(sFig.sObject(intObject).handles.lines)
					delete(sFig.sObject(intObject).handles.lines(p));
				end
			else
				p = [];
			end
			if isempty(p) && isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX)
				delete(sFig.sObject(intObject).handles.marker);
			end
			if isfield(sFig.sObject(intObject).handles,'text') && ~isempty(sFig.sObject(intObject).handles.text)
				delete(sFig.sObject(intObject).handles.text);
				sFig.sObject(intObject).handles.text = [];
			end
			sFig.sObject(intObject).drawn = 0;

			%remove mask
			sDC.ROI(intObject).matMask = [];
		end
		
		%keep all others
		sDC.ROI = sDC.ROI(~vecRemIndex);
		sFig.sObject = sFig.sObject(~vecRemIndex);
	end
	%set msg
	cellText{1} = ['Removed ' num2str(length(vecRemObjects)) ' objects'] ;
	DC_updateTextInformation(cellText)
end
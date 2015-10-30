function intThisType = DC_getTotalType()
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
	%get structures
	global sRec;
	global sDC;
	global sFig;
	
	%check if data has been loaded
	if isempty(sDC) || isempty(sFig) || isempty(sRec)
		return;
	else
		try
			%get current image
			intImSelected = get(sFig.ptrListSelectImage,'Value');
		catch %#ok<CTCH>
			return;
		end
	end
	
	%get current type
	intType = get(sFig.ptrListCellType, 'Value');
	
	%get total number for this type
	intObjects = numel(sDC.ROI);
	intThisType = 0;
	for intObject=1:intObjects
		if intType == sDC.ROI(intObject).intType
			intThisType = intThisType + 1;
		end
	end
end


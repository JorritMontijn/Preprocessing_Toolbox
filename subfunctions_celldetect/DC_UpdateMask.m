function DC_UpdateMask(intObject,matMask)
	%removes objects from sDC structure
	%DC_UpdateMask(intObject,matMask)
	
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
	
	%perform checks
	if sum(matMask(:)) == 0
		if ~isnan(sDC.ROI(intObject).intCenterY) && ~isnan(sDC.ROI(intObject).intCenterX)
			%set region around centroid as mask
			sDC.ROI(intObject).matMask = matMask;
			sDC.ROI(intObject).matMask([-1:1]+round(sDC.ROI(intObject).intCenterY),[-1:1]+round(sDC.ROI(intObject).intCenterX)) = true;
		else
			%otherwise, don't update
			warning([mfilename ':NoCentroid'],sprintf('Empty mask and empty centroid for object %d; skipping...',intObject));
		end
	else
		%set mask
		sDC.ROI(intObject).matMask = matMask;
	end
	%set region around centroid as mask
	[vecYX] = calcCenterOfMass(sDC.ROI(intObject).matMask);
	sDC.ROI(intObject).intCenterX = vecYX(2);
	sDC.ROI(intObject).intCenterY = vecYX(1);		
end
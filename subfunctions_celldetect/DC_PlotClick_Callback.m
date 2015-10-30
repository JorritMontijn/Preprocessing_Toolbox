function DC_PlotClick_Callback(hObject, eventdata) %#ok<INUSD>
	
	%get globals
	global sFig;
	global sDC;
	
	%get type of click
	strClickType = get(gcf, 'SelectionType');
	if strcmpi(strClickType,'alt') && ~isempty(sFig.vecSelectedObjects)
		%lock GUI
		DC_lock(sFig);
		
		%remove all selected objects
		DC_unselectAll();
		
		%set selection & update text
		DC_updateTextInformation;
		
		%redraw
		DC_redraw(0);
		
		%unlock GUI
		DC_unlock(sFig);
	elseif strcmpi(strClickType,'extend')
		%lock GUI
		DC_lock(sFig);
		
		%get click area
		rect = getrect;
		intMinX=rect(1);
		intMinY=rect(2);
		intMaxX=rect(1) + rect(3);
		intMaxY=rect(2) + rect(4);
		
		% Get margins of screen
		axesSize = get(sFig.ptrAxesHandle, 'Position');
		xZero = axesSize(1);
		yZero = axesSize(2);
		width = axesSize(3);
		height = axesSize(4);
		
		% get margins of plot-window
		sizeI = size(sFig.imCurrent);
		screenWidth = sizeI(2);
		screenHeight = sizeI(1);
		
		% Calculate the coordinates of the click related to the image
		% coordinates
		dblMagnification = str2double(get(sFig.ptrEditMagnification,'String'))/100;
		dblRealMinX = ( ((intMinX - xZero) * screenWidth ) / width  )/dblMagnification;
		dblRealMinY = ( ((intMinY - yZero) * screenHeight) / height )/dblMagnification;
		dblRealMaxX = ( ((intMaxX - xZero) * screenWidth ) / width  )/dblMagnification;
		dblRealMaxY = ( ((intMaxY - yZero) * screenHeight) / height )/dblMagnification;
		intRealMinX = max(round(dblRealMinX),1);
		intRealMinY = max(round(dblRealMinY),1);
		intRealMaxX = min(round(dblRealMaxX),floor(screenWidth/dblMagnification));
		intRealMaxY = min(round(dblRealMaxY),floor(screenHeight/dblMagnification));
		
		%check draw type
		intDrawType = get(sFig.ptrListSelectDrawType,'Value');
		cellDrawTypes = get(sFig.ptrListSelectDrawType,'String');
		strDrawType = cellDrawTypes{intDrawType};
		boolDrawBorder = strcmpi(strDrawType,'Border');
		
		% check if cell was clicked
		boolChange = false;
		intObjects = numel(sDC.ROI);
		for intObject = 1:intObjects
			if isfield(sDC.ROI,'matMask') && ~isempty(sDC.ROI(intObject).matMask) && sum(sum(sDC.ROI(intObject).matMask(intRealMinY:intRealMaxY,intRealMinX:intRealMaxX))) > 0
				
				%set select flag
				vecSelectIndex = sFig.vecSelectedObjects == intObject;
				if max(vecSelectIndex)
					sFig.vecSelectedObjects = sFig.vecSelectedObjects(~vecSelectIndex);
				else
					sFig.vecSelectedObjects(end+1) = intObject;
				end
				
				% remove old drawing and set drawn flag
				if boolDrawBorder
					for p = 1:length(sFig.sObject(intObject).handles.lines)
						delete(sFig.sObject(intObject).handles.lines(p));
					end
				else
					p = [];
				end
				if isempty(p) && isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX)
					delete(sFig.sObject(intObject).handles.marker);
				end
				sFig.sObject(intObject).drawn = 0;
				
				% end loop
				boolChange = true;
			end
		end
		
		%check if fig needs updating
		if boolChange
			%update information window
			DC_updateTextInformation;
			
			%redraw outlines
			DC_redraw(0);
		end
		
		%unlock GUI
		DC_unlock(sFig);
	elseif strcmpi(strClickType,'normal')
		% Get margins of screen
		axesSize = get(sFig.ptrAxesHandle, 'Position');
		xZero = axesSize(1);
		yZero = axesSize(2);
		width = axesSize(3);
		height = axesSize(4);
		
		% get margins of plot-window
		sizeI = size(sFig.imCurrent);
		screenWidth = sizeI(2);
		screenHeight = sizeI(1);
		
		% get location of the click
		XY = get(hObject, 'CurrentPoint');
		X = XY(1);
		Y = XY(2);
		
		% Calculate the coordinates of the click related to the image
		% coordinates
		dblMagnification = str2double(get(sFig.ptrEditMagnification,'String'))/100;
		dblRealX = ( ((X - xZero) * screenWidth ) / width  )/dblMagnification;
		dblRealMinY = ( ((Y - yZero) * screenHeight) / height )/dblMagnification;
		intRealMinX = round(dblRealX);
		intRealMinY = round(dblRealMinY);
		
		%check draw type
		intDrawType = get(sFig.ptrListSelectDrawType,'Value');
		cellDrawTypes = get(sFig.ptrListSelectDrawType,'String');
		strDrawType = cellDrawTypes{intDrawType};
		boolDrawBorder = strcmpi(strDrawType,'Border');
		
		% check if cell was clicked
		boolChange = false;
		intObjects = numel(sDC.ROI);
		for intObject = 1:intObjects
			if isfield(sDC.ROI,'matMask') && ~isempty(sDC.ROI(intObject).matMask) && sDC.ROI(intObject).matMask(intRealMinY,intRealMinX)
				
				%set select flag
				vecSelectIndex = sFig.vecSelectedObjects == intObject;
				if max(vecSelectIndex)
					sFig.vecSelectedObjects = sFig.vecSelectedObjects(~vecSelectIndex);
				else
					sFig.vecSelectedObjects(end+1) = intObject;
				end
				
				% remove old drawing and set drawn flag
				if boolDrawBorder
					for p = 1:length(sFig.sObject(intObject).handles.lines)
						delete(sFig.sObject(intObject).handles.lines(p));
					end
				else
					p = [];
				end
				if isempty(p) && isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX)
					delete(sFig.sObject(intObject).handles.marker);
				end
				sFig.sObject(intObject).drawn = 0;
				
				% end loop
				boolChange = true;
				break;
			end
		end
		
		%check if fig needs updating
		if boolChange
			%update information window
			DC_updateTextInformation;
			
			%redraw outlines
			DC_redraw(0);
		end
	end
end
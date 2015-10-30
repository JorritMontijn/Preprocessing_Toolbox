function DC_moveROIs
	%globals
	global sDC;
	global sFig;
	global sMoveROI;
	
	%check if data has been loaded
	if isempty(sDC) || isempty(sFig)
		return;
	else
		try
			intType = get(sFig.ptrListCellType, 'Value');
		catch
			return;
		end
	end
	
	%get global data
	sMoveROI.intROIDisplacementX = 0;
	sMoveROI.intROIDisplacementY = 0;
	sMoveROI.boolUpdated = false;
	sMoveROI.boolRunning = true;
	
	%define variables
	intROIDisplacementX = 0;
	intROIDisplacementY = 0;
	intROIAssignedDispX = 0;
	intROIAssignedDispY = 0;
	
	%put GUI on top
	figure(sFig.ptrMainGUI);
	figure(sFig.ptrWindowHandle);
	
	%run GUI
	DC_MoveROIsGUI;
	
	%get displacement value
	while sMoveROI.boolRunning
		% if updated
		if sMoveROI.boolUpdated
			%unset update flag
			sMoveROI.boolUpdated = false;
			
			%get updated positions
			intROIDisplacementX = sMoveROI.intROIDisplacementX;
			intROIDisplacementY = sMoveROI.intROIDisplacementY;
			
			cellText = {'Computing...'};
			DC_updateTextInformation(cellText);
		
			%get incrementally required displacement
			intMoveX = intROIDisplacementX - intROIAssignedDispX;
			intMoveY = intROIDisplacementY - intROIAssignedDispY;
			intROIAssignedDispX = intROIDisplacementX;
			intROIAssignedDispY = intROIDisplacementY;
			
			%assign displacement to all objects
			intObjects = numel(sDC.ROI);
			for intObject = 1:intObjects
				if ismember(intObject,sFig.vecSelectedObjects) || isempty(sFig.vecSelectedObjects)
					if ~isempty(sDC.ROI(intObject).matMask)
						%translocate perimeter
						%sDC.ROI(intObject).matPerimeter(:, 1) = sDC.ROI(intObject).matPerimeter(:, 1) + intMoveX;
						%sDC.ROI(intObject).matPerimeter(:, 2) = sDC.ROI(intObject).matPerimeter(:, 2) + intMoveY;
						
						sDC.ROI(intObject).matMask = circshift(sDC.ROI(intObject).matMask,[intMoveY intMoveX]);
					end
					if isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX)
						%draw center
						sDC.ROI(intObject).intCenterX = sDC.ROI(intObject).intCenterX + intMoveX;
						sDC.ROI(intObject).intCenterY = sDC.ROI(intObject).intCenterY + intMoveY;
					end
					
					% remove old drawing and set drawn flag
					if isfield(sFig.sObject(intObject).handles,'text') && ~isempty(sFig.sObject(intObject).handles.text)
						delete(sFig.sObject(intObject).handles.text);
						sFig.sObject(intObject).handles.text = [];
					end
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
					sFig.sObject(intObject).drawn = 0;
				end
			end
			
			%redraw
			DC_redraw(0);
			
			%set msg
			cellText{1} = ['Help info ROI movement'] ;
			cellText{2} = [''] ;
			cellText{3} = sprintf('Current displacement; x: %.0f; y=%.0f',sMoveROI.intROIDisplacementX,sMoveROI.intROIDisplacementY) ;
			cellText{4} = [''] ;
			cellText{5} = ['Will move all objects if none are selected; or will move only selected objects'] ;
			DC_updateTextInformation(cellText)
		else
			pause(0.1);
		end
	end
	%get incrementally required displacement
	intMoveX = intROIDisplacementX - intROIAssignedDispX;
	intMoveY = intROIDisplacementY - intROIAssignedDispY;
	intROIAssignedDispX = intROIDisplacementX;
	intROIAssignedDispY = intROIDisplacementY;
	
	if intMoveX ~= 0 || intMoveY ~= 0
		%assign displacement to all objects
		intObjects = numel(sDC.ROI);
		for intObject = 1:intObjects
			if ismember(intObject,sFig.vecSelectedObjects) || isempty(sFig.vecSelectedObjects)
				if ~isempty(sDC.ROI(intObject).matMask)
					%translocate perimeter
					%sDC.ROI(intObject).matPerimeter(:, 1) = sDC.ROI(intObject).matPerimeter(:, 1) + intMoveX;
					%sDC.ROI(intObject).matPerimeter(:, 2) = sDC.ROI(intObject).matPerimeter(:, 2) + intMoveY;
						
					sDC.ROI(intObject).matMask = circshift(sDC.ROI(intObject).matMask,[intMoveY intMoveX]);
				end
				if isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX)
					%draw center
					sDC.ROI(intObject).intCenterX = sDC.ROI(intObject).intCenterX + intMoveX;
					sDC.ROI(intObject).intCenterY = sDC.ROI(intObject).intCenterY + intMoveY;
				end
				
				% remove old drawing and set drawn flag
				for p = 1:length(sFig.sObject(intObject).handles.lines)
					delete(sFig.sObject(intObject).handles.lines(p));
				end
				if isempty(p) && isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX)
					delete(sFig.sObject(intObject).handles.marker);
				end
				sFig.sObject(intObject).drawn = 0;
			end
		end
		
		%redraw
		DC_redraw(0);
	end
	
	if isempty(sFig.vecSelectedObjects)
		%save to structure
		sDC.metaData.intROIDisplacementX = intROIDisplacementX;
		sDC.metaData.intROIDisplacementY = intROIDisplacementY;
		sDC.metaData.intROIAssignedDispX = intROIAssignedDispX;
		sDC.metaData.intROIAssignedDispY = intROIAssignedDispY;
	end
	
	%set msg
	clear cellText;
	cellText{1} = ['Completed ROI location movement'] ;
	cellText{2} = [''] ;
	cellText{3} = sprintf('Current displacement; x: %.0f; y=%.0f',intROIDisplacementX,intROIDisplacementY) ;
	DC_updateTextInformation(cellText)
end